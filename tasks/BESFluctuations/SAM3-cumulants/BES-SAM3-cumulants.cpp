#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdio>
#include <cassert>

#include "HRGBase.h"
#include "HRGEV.h"
#include "HRGFit.h"
#include "HRGVDW.h"
#include "HRGEventGenerator.h"

#include "ThermalFISTConfig.h"

#include "HypersurfaceReader.h"
#include "FistSamplerParameters.h"
#include "FistSamplerHelperFunctions.h"

#include "FistSamplerConfig.h"

#include "sample-moments/include/NumberStatistics.h"
#include "sample-moments/include/TwoNumberStatistics.h"

using namespace std;

#ifdef ThermalFIST_USENAMESPACE
using namespace thermalfist;
#endif

using namespace FistSampler;
using namespace SampleMoments;


// ============================================================
// Event processor: collects GCE cumulants needed for SAM-3.0
//
// Observables (d=2):  N_p, N_pbar  in STAR acceptance
//                     (|y| < y_cut, pTmin < pT < pTmax)
// Conserved charges (s=3):  B_tot, Q_tot, S_tot  (4pi)
//
// Computes all 2nd-order (cross-)cumulants as functions of y_cut
// ============================================================
class EventsProcessorSAM3 {
public:
  int m_nsubs;
  double m_dY;
  double m_pTmin, m_pTmax;
  long long nevents;

  // --- Acceptance-dependent (per rapidity bin) ---

  // p/pbar joint statistics:  var(Np), var(Npbar), cov(Np,Npbar)
  vector<TwoNumberStatistics> statsPPbar;

  // Cross-cumulants of accepted p with 4pi charges
  vector<TwoNumberStatistics> statsPB, statsPQ, statsPS;

  // Cross-cumulants of accepted pbar with 4pi charges
  vector<TwoNumberStatistics> statsPbarB, statsPbarQ, statsPbarS;

  // --- 4pi charge statistics (acceptance-independent) ---

  // Charge pairs: cov(B,Q), cov(B,S), cov(Q,S)
  TwoNumberStatistics statsBQ, statsBS, statsQS;

  // Individual charges: var(B), var(Q), var(S)
  NumberStatistics statsB, statsQ, statsS;


  EventsProcessorSAM3(int nsubs, double dY, double pTmin, double pTmax)
    : m_nsubs(nsubs), m_dY(dY), m_pTmin(pTmin), m_pTmax(pTmax), nevents(0)
  {
    int nbins = m_nsubs / 2;
    statsPPbar.resize(nbins);
    statsPB.resize(nbins);
    statsPQ.resize(nbins);
    statsPS.resize(nbins);
    statsPbarB.resize(nbins);
    statsPbarQ.resize(nbins);
    statsPbarS.resize(nbins);
  }


  /// Process a primordial event: compute 4pi charges before decays,
  /// then perform decays and count accepted p/pbar from final state.
  /// This ensures conserved charges (especially S) are not affected by
  /// weak decays that can violate strangeness conservation.
  void ProcessEvent(const SimpleEvent& evt, ThermalParticleSystem& TPS, int decays) {
    nevents++;

    int nbins = m_nsubs / 2;
    vector<int> Np(nbins, 0), Npbar(nbins, 0);
    int Btot = 0, Qtot = 0, Stot = 0;

    // 4pi conserved charges from primordial particles (before decays)
    for (const SimpleParticle& part : evt.Particles) {
      const ThermalParticle& prop = TPS.ParticleByPDG(part.PDGID);
      Btot += prop.BaryonCharge();
      Qtot += prop.ElectricCharge();
      Stot += prop.Strangeness();
    }

    // Perform decays
    SimpleEvent evt_final = evt;
    if (decays != 0) {
      evt_final = EventGeneratorBase::PerformDecays(evt, &TPS);
    }

    // Accepted protons/antiprotons from final-state particles (after decays)
    for (const SimpleParticle& part : evt_final.Particles) {
      double pT = part.GetPt();
      if (pT >= m_pTmin && pT <= m_pTmax) {
        int tindY = floor(abs(part.GetY()) / m_dY);
        if (tindY >= 0 && tindY < nbins) {
          if (part.PDGID == 2212)
            Np[tindY]++;
          if (part.PDGID == -2212)
            Npbar[tindY]++;
        }
      }
    }

    // Prefix sums: bin i -> |y| < (i+1)*dY
    for (int i = 1; i < nbins; ++i) {
      Np[i] += Np[i - 1];
      Npbar[i] += Npbar[i - 1];
    }

    // Fill acceptance-dependent statistics
    for (int isub = 0; isub < nbins; ++isub) {
      statsPPbar[isub].AddObservation(Np[isub], Npbar[isub]);

      statsPB[isub].AddObservation(Np[isub], Btot);
      statsPQ[isub].AddObservation(Np[isub], Qtot);
      statsPS[isub].AddObservation(Np[isub], Stot);

      statsPbarB[isub].AddObservation(Npbar[isub], Btot);
      statsPbarQ[isub].AddObservation(Npbar[isub], Qtot);
      statsPbarS[isub].AddObservation(Npbar[isub], Stot);
    }

    // 4pi charge correlations
    statsBQ.AddObservation(Btot, Qtot);
    statsBS.AddObservation(Btot, Stot);
    statsQS.AddObservation(Qtot, Stot);

    statsB.AddObservation(Btot);
    statsQ.AddObservation(Qtot);
    statsS.AddObservation(Stot);
  }
};


// ============================================================
// File writer
// ============================================================
void WriteToFile(const string& prefix, EventsProcessorSAM3& stats) {
  ofstream fout;
  int w = 15;

  fout.open(prefix + ".SAM3-cumulants.dat");

  // ---- Header comments ----
  fout << "# GCE cumulants for SAM-3.0 CE corrections (Eq. 53)" << endl;
  fout << "# kappa^ce_{ij}(X) = kappa^gce_{ij;0} - kappa^gce_{i;a} (K^{-1})^{ab} kappa^gce_{j;b}" << endl;
  fout << "# Observables X = (N_p, N_pbar),  Conserved charges = (B, Q, S) in 4pi" << endl;
  fout << "# Events: " << stats.nevents << endl;
  fout << "# pT cuts: " << stats.m_pTmin << " < pT < " << stats.m_pTmax << " GeV/c" << endl;
  fout << "# 4pi charges computed before decays" << endl;
  fout << "#" << endl;

  // ---- Compute 4pi charge quantities (acceptance-independent) ----
  double meanB  = stats.statsB.GetMean();
  double meanQ  = stats.statsQ.GetMean();
  double meanS  = stats.statsS.GetMean();

  double varB  = stats.statsB.GetCentralMoment(2);
  double varQ  = stats.statsQ.GetCentralMoment(2);
  double varS  = stats.statsS.GetCentralMoment(2);
  double covBQ = stats.statsBQ.GetJointCentralMoment(1, 1);
  double covBS = stats.statsBS.GetJointCentralMoment(1, 1);
  double covQS = stats.statsQS.GetJointCentralMoment(1, 1);

  double meanBe = stats.statsB.GetMeanError();
  double meanQe = stats.statsQ.GetMeanError();
  double meanSe = stats.statsS.GetMeanError();
  double varBe  = stats.statsB.GetCentralMomentError(2);
  double varQe  = stats.statsQ.GetCentralMomentError(2);
  double varSe  = stats.statsS.GetCentralMomentError(2);
  double covBQe = stats.statsBQ.GetJointCentralMomentError(1, 1);
  double covBSe = stats.statsBS.GetJointCentralMomentError(1, 1);
  double covQSe = stats.statsQS.GetJointCentralMomentError(1, 1);

  // ---- Print K matrix in header for reference ----
  fout << "# 4pi charge means:  <B> = " << meanB << " +/- " << meanBe
       << "   <Q> = " << meanQ << " +/- " << meanQe
       << "   <S> = " << meanS << " +/- " << meanSe << endl;
  fout << "#" << endl;
  fout << "# 4pi charge covariance matrix K_ab (value +/- error):" << endl;
  fout << "#       B                Q                S" << endl;
  fout << "# B  " << setw(w) << varB  << " +/- " << setw(w) << varBe
       << "  " << setw(w) << covBQ << " +/- " << setw(w) << covBQe
       << "  " << setw(w) << covBS << " +/- " << setw(w) << covBSe << endl;
  fout << "# Q  " << setw(w) << covBQ << " +/- " << setw(w) << covBQe
       << "  " << setw(w) << varQ  << " +/- " << setw(w) << varQe
       << "  " << setw(w) << covQS << " +/- " << setw(w) << covQSe << endl;
  fout << "# S  " << setw(w) << covBS << " +/- " << setw(w) << covBSe
       << "  " << setw(w) << covQS << " +/- " << setw(w) << covQSe
       << "  " << setw(w) << varS  << " +/- " << setw(w) << varSe  << endl;
  fout << "#" << endl;

  // ---- Column header line ----
  fout << setw(w) << "ycut"
       << setw(w) << "<Np>"         << setw(w) << "<Np>_err"
       << setw(w) << "<Npbar>"      << setw(w) << "<Npbar>_err"
       << setw(w) << "var(Np)"      << setw(w) << "var(Np)_err"
       << setw(w) << "var(Npbar)"   << setw(w) << "var(Npb)_err"
       << setw(w) << "cov(Np,Npb)"  << setw(w) << "cov(NpNpb)er"
       << setw(w) << "cov(Np,B)"    << setw(w) << "cov(NpB)_err"
       << setw(w) << "cov(Np,Q)"    << setw(w) << "cov(NpQ)_err"
       << setw(w) << "cov(Np,S)"    << setw(w) << "cov(NpS)_err"
       << setw(w) << "cov(Npb,B)"   << setw(w) << "cov(NpbB)err"
       << setw(w) << "cov(Npb,Q)"   << setw(w) << "cov(NpbQ)err"
       << setw(w) << "cov(Npb,S)"   << setw(w) << "cov(NpbS)err"
       << setw(w) << "<B>"          << setw(w) << "<Q>"          << setw(w) << "<S>"
       << setw(w) << "var(B)"       << setw(w) << "cov(B,Q)"     << setw(w) << "cov(B,S)"
       << setw(w) << "var(Q)"       << setw(w) << "cov(Q,S)"     << setw(w) << "var(S)"
       << endl;

  // ---- Data rows ----
  int nbins = stats.m_nsubs / 2;
  for (int isub = 0; isub < nbins; ++isub) {
    double ycut = (isub + 1) * stats.m_dY;

    auto& ppbar = stats.statsPPbar[isub];

    double meanNp     = ppbar.GetMean1();
    double meanNpErr  = ppbar.GetMean1Error();
    double meanNpb    = ppbar.GetMean2();
    double meanNpbErr = ppbar.GetMean2Error();

    double varNp        = ppbar.GetJointCumulant(2, 0);
    double varNpErr     = ppbar.GetJointCumulantError(2, 0);
    double varNpb       = ppbar.GetJointCumulant(0, 2);
    double varNpbErr    = ppbar.GetJointCumulantError(0, 2);
    double covNpNpb     = ppbar.GetJointCumulant(1, 1);
    double covNpNpbErr  = ppbar.GetJointCumulantError(1, 1);

    double covNpB       = stats.statsPB[isub].GetJointCumulant(1, 1);
    double covNpBErr    = stats.statsPB[isub].GetJointCumulantError(1, 1);
    double covNpQ       = stats.statsPQ[isub].GetJointCumulant(1, 1);
    double covNpQErr    = stats.statsPQ[isub].GetJointCumulantError(1, 1);
    double covNpS       = stats.statsPS[isub].GetJointCumulant(1, 1);
    double covNpSErr    = stats.statsPS[isub].GetJointCumulantError(1, 1);

    double covNpbB      = stats.statsPbarB[isub].GetJointCumulant(1, 1);
    double covNpbBErr   = stats.statsPbarB[isub].GetJointCumulantError(1, 1);
    double covNpbQ      = stats.statsPbarQ[isub].GetJointCumulant(1, 1);
    double covNpbQErr   = stats.statsPbarQ[isub].GetJointCumulantError(1, 1);
    double covNpbS      = stats.statsPbarS[isub].GetJointCumulant(1, 1);
    double covNpbSErr   = stats.statsPbarS[isub].GetJointCumulantError(1, 1);

    fout << setw(w) << ycut
         << setw(w) << meanNp     << setw(w) << meanNpErr
         << setw(w) << meanNpb    << setw(w) << meanNpbErr
         << setw(w) << varNp      << setw(w) << varNpErr
         << setw(w) << varNpb     << setw(w) << varNpbErr
         << setw(w) << covNpNpb   << setw(w) << covNpNpbErr
         << setw(w) << covNpB     << setw(w) << covNpBErr
         << setw(w) << covNpQ     << setw(w) << covNpQErr
         << setw(w) << covNpS     << setw(w) << covNpSErr
         << setw(w) << covNpbB    << setw(w) << covNpbBErr
         << setw(w) << covNpbQ    << setw(w) << covNpbQErr
         << setw(w) << covNpbS    << setw(w) << covNpbSErr
         << setw(w) << meanB      << setw(w) << meanQ      << setw(w) << meanS
         << setw(w) << varB       << setw(w) << covBQ      << setw(w) << covBS
         << setw(w) << varQ       << setw(w) << covQS      << setw(w) << varS
         << endl;
  }

  fout.close();
}


// ============================================================
// Write SAM-3.0 CE-corrected 2nd-order cumulants
//
// Applies Eq. (53):
//   kappa^ce_{ij} = kappa^gce_{ij} - kappa^gce_{i,a} (K^{-1})^{ab} kappa^gce_{j,b}
//
// for different canonical scenarios:
//   B-can, Q-can, S-can, BQ-can, BS-can, BQS-can
// ============================================================
void WriteSAM3CorrectedFile(const string& prefix, EventsProcessorSAM3& stats) {
  ofstream fout;
  int w = 15;

  fout.open(prefix + ".SAM3-corrected.dat");

  fout << "# SAM-3.0 CE-corrected 2nd-order cumulants for p, pbar" << endl;
  fout << "# kappa^ce_{ij} = kappa^gce_{ij} - kappa^gce_{i,a} (K^{-1})^{ab} kappa^gce_{j,b}" << endl;
  fout << "# Observables: X = (N_p, N_pbar) in STAR acceptance" << endl;
  fout << "# Conserved charges: B, Q, S  (4pi, before decays)" << endl;
  fout << "# Events: " << stats.nevents << endl;
  fout << "# pT cuts: " << stats.m_pTmin << " < pT < " << stats.m_pTmax << " GeV/c" << endl;
  fout << "#" << endl;

  // ---- 4pi charge covariance matrix K_{ab} ----
  double varB  = stats.statsB.GetCentralMoment(2);
  double varQ  = stats.statsQ.GetCentralMoment(2);
  double varS  = stats.statsS.GetCentralMoment(2);
  double covBQ = stats.statsBQ.GetJointCentralMoment(1, 1);
  double covBS = stats.statsBS.GetJointCentralMoment(1, 1);
  double covQS = stats.statsQS.GetJointCentralMoment(1, 1);

  // Full 3x3 K matrix: indices 0=B, 1=Q, 2=S
  double K[3][3] = {
    {varB,  covBQ, covBS},
    {covBQ, varQ,  covQS},
    {covBS, covQS, varS}
  };

  fout << "# K matrix (4pi charge covariances):" << endl;
  fout << "#       B              Q              S" << endl;
  fout << "# B  " << setw(w) << varB  << setw(w) << covBQ << setw(w) << covBS << endl;
  fout << "# Q  " << setw(w) << covBQ << setw(w) << varQ  << setw(w) << covQS << endl;
  fout << "# S  " << setw(w) << covBS << setw(w) << covQS << setw(w) << varS  << endl;
  fout << "#" << endl;

  // ---- Helper: invert symmetric NxN submatrix (N=1,2,3) ----
  // idx: which charge indices from {0,1,2} are conserved
  // Returns K^{-1} as flat vector of size N*N (row-major)
  auto invertSubmatrix = [&K](const vector<int>& idx) -> vector<double> {
    int N = (int)idx.size();
    if (N == 1) {
      return {1.0 / K[idx[0]][idx[0]]};
    }
    else if (N == 2) {
      double a = K[idx[0]][idx[0]], b = K[idx[0]][idx[1]];
      double d = K[idx[1]][idx[1]];
      double det = a * d - b * b;
      return {d / det, -b / det, -b / det, a / det};
    }
    else {
      // 3x3 inversion via cofactors
      double M[3][3];
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          M[i][j] = K[idx[i]][idx[j]];

      double det = M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
                 - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
                 + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]);

      vector<double> inv(9);
      inv[0] =  (M[1][1]*M[2][2] - M[1][2]*M[2][1]) / det;
      inv[1] = -(M[0][1]*M[2][2] - M[0][2]*M[2][1]) / det;
      inv[2] =  (M[0][1]*M[1][2] - M[0][2]*M[1][1]) / det;
      inv[3] = -(M[1][0]*M[2][2] - M[1][2]*M[2][0]) / det;
      inv[4] =  (M[0][0]*M[2][2] - M[0][2]*M[2][0]) / det;
      inv[5] = -(M[0][0]*M[1][2] - M[0][2]*M[1][0]) / det;
      inv[6] =  (M[1][0]*M[2][1] - M[1][1]*M[2][0]) / det;
      inv[7] = -(M[0][0]*M[2][1] - M[0][1]*M[2][0]) / det;
      inv[8] =  (M[0][0]*M[1][1] - M[0][1]*M[1][0]) / det;
      return inv;
    }
  };

  // ---- Compute correction: Delta_{ij} = sum_{a,b} cp_i[a] * Kinv[a,b] * cp_j[b] ----
  auto computeCorrection = [](const vector<double>& Kinv, const vector<int>& idx,
                               const double cp_i[3], const double cp_j[3]) -> double {
    int N = (int)idx.size();
    double result = 0.;
    for (int a = 0; a < N; ++a)
      for (int b = 0; b < N; ++b)
        result += cp_i[idx[a]] * Kinv[a * N + b] * cp_j[idx[b]];
    return result;
  };

  // ---- Define canonical scenarios ----
  struct Scenario {
    string name;
    vector<int> idx;  // which charges are conserved: 0=B, 1=Q, 2=S
  };

  vector<Scenario> scenarios = {
    {"Bcan",   {0}},
    {"Qcan",   {1}},
    {"Scan",   {2}},
    {"BQcan",  {0, 1}},
    {"BScan",  {0, 2}},
    {"BQScan", {0, 1, 2}}
  };

  // Precompute K^{-1} for each scenario
  vector<vector<double>> Kinvs;
  for (auto& sc : scenarios)
    Kinvs.push_back(invertSubmatrix(sc.idx));

  // Print K^{-1} matrices in header
  for (int isc = 0; isc < (int)scenarios.size(); ++isc) {
    auto& sc = scenarios[isc];
    auto& Kinv = Kinvs[isc];
    int N = (int)sc.idx.size();
    const char* chargeLabel[] = {"B", "Q", "S"};
    fout << "# K^{-1} for " << sc.name << ":";
    for (int a = 0; a < N; ++a)
      for (int b = 0; b < N; ++b)
        fout << "  " << chargeLabel[sc.idx[a]] << chargeLabel[sc.idx[b]] << "=" << Kinv[a * N + b];
    fout << endl;
  }
  fout << "#" << endl;

  // ---- Column header line ----
  fout << setw(w) << "ycut"
       << setw(w) << "<Np>"
       << setw(w) << "<Npbar>"
       << setw(w) << "k2p_gce"
       << setw(w) << "k2pbar_gce"
       << setw(w) << "k11ppb_gce";

  for (auto& sc : scenarios) {
    fout << setw(w) << ("k2p_" + sc.name)
         << setw(w) << ("k2pb_" + sc.name)
         << setw(w) << ("k11_" + sc.name);
  }
  fout << endl;

  // ---- Data rows ----
  int nbins = stats.m_nsubs / 2;
  for (int isub = 0; isub < nbins; ++isub) {
    double ycut = (isub + 1) * stats.m_dY;

    auto& ppbar = stats.statsPPbar[isub];

    double meanNp   = ppbar.GetMean1();
    double meanNpb  = ppbar.GetMean2();

    // GCE 2nd-order cumulants
    double varNp    = ppbar.GetJointCumulant(2, 0);
    double varNpb   = ppbar.GetJointCumulant(0, 2);
    double covNpNpb = ppbar.GetJointCumulant(1, 1);

    // Cross-cumulants with charges: index 0=B, 1=Q, 2=S
    double cpNp[3] = {
      stats.statsPB[isub].GetJointCumulant(1, 1),
      stats.statsPQ[isub].GetJointCumulant(1, 1),
      stats.statsPS[isub].GetJointCumulant(1, 1)
    };
    double cpNpb[3] = {
      stats.statsPbarB[isub].GetJointCumulant(1, 1),
      stats.statsPbarQ[isub].GetJointCumulant(1, 1),
      stats.statsPbarS[isub].GetJointCumulant(1, 1)
    };

    fout << setw(w) << ycut
         << setw(w) << meanNp
         << setw(w) << meanNpb
         << setw(w) << varNp
         << setw(w) << varNpb
         << setw(w) << covNpNpb;

    for (int isc = 0; isc < (int)scenarios.size(); ++isc) {
      auto& sc = scenarios[isc];
      auto& Kinv = Kinvs[isc];

      // CE-corrected = GCE - correction
      double dVarNp    = computeCorrection(Kinv, sc.idx, cpNp,  cpNp);
      double dVarNpb   = computeCorrection(Kinv, sc.idx, cpNpb, cpNpb);
      double dCovNpNpb = computeCorrection(Kinv, sc.idx, cpNp,  cpNpb);

      fout << setw(w) << (varNp    - dVarNp)
           << setw(w) << (varNpb   - dVarNpb)
           << setw(w) << (covNpNpb - dCovNpNpb);
    }
    fout << endl;
  }

  fout.close();
}


// ============================================================
// Main
// ============================================================
FistSamplerParameters run_parameters;

int nsubs = 200;       // Number of subintervals for Y
double dY = 0.05;     // Rapidity cut step
double pTmin = 0.4;   // STAR proton pT min [GeV/c]
double pTmax = 2.0;   // STAR proton pT max [GeV/c]


int main(int argc, char* argv[]) {

  cout << "Running FIST sampler version " << FistSampler_VERSION_MAJOR << "." << FistSampler_VERSION_MINOR << endl;
  cout << "SAM-3.0 cumulants calculation" << endl << endl;

  run_parameters.parameters["ecm"] = 7.7;

  run_parameters.ReadParametersFromCommandLine(argc, argv);

  // nsubs
  if (run_parameters.parameters.count("nsubs") > 0)
    nsubs = lround(run_parameters.parameters["nsubs"]);

  // pT cuts
  if (run_parameters.parameters.count("pTmin") > 0)
    pTmin = run_parameters.parameters["pTmin"];
  if (run_parameters.parameters.count("pTmax") > 0)
    pTmax = run_parameters.parameters["pTmax"];

  // Prefix with ensemble suffix (e.g. GCE, B, BS, BQS)
  string prefix = run_parameters.output_file;
  size_t lastindex = prefix.find_last_of(".");
  prefix = prefix.substr(0, lastindex);

  {
    int Bcan = lround(run_parameters.parameters["Bcanonical"]);
    int Qcan = lround(run_parameters.parameters["Qcanonical"]);
    int Scan = lround(run_parameters.parameters["Scanonical"]);
    string ensemble_suffix;
    if (!Bcan && !Qcan && !Scan)
      ensemble_suffix = ".GCE";
    else {
      ensemble_suffix = ".";
      if (Bcan) ensemble_suffix += "B";
      if (Qcan) ensemble_suffix += "Q";
      if (Scan) ensemble_suffix += "S";
    }
    prefix += ensemble_suffix;
  }

  // Output the values of all the parameters used
  run_parameters.OutputParameters();

  // Set the random seed
  RandomGenerators::SetSeed(run_parameters.randomseed);

  int fist_sampler_mode = lround(run_parameters.parameters["fist_sampler_mode"]);
  if (fist_sampler_mode < 0 || fist_sampler_mode > 2) {
    cout << "fist_sampler_mode of " << fist_sampler_mode << " is unsupported! Aborting..." << endl;
    exit(1);
  }

  // Cooper-Frye hypersurface
  ParticlizationHypersurface hypersurface;

  if (fist_sampler_mode == 0) {
    ReadHypersurfaceFromFile(run_parameters, hypersurface);
    if (hypersurface.size() == 0) {
      cout << "Empty hypersurface! Aborting..." << endl;
      exit(1);
    }
  }
  else if (fist_sampler_mode == 1) {
    CreateSiemensRasmussenHubbleHypersurface(run_parameters, hypersurface);
  }

  cout << "Initializing event generator..." << endl;
  EventGeneratorBase* evtgen;

  if (fist_sampler_mode == 0 || fist_sampler_mode == 1) {
    evtgen = CreateEventGeneratorFromHypersurface(run_parameters, hypersurface);
  }
  else if (fist_sampler_mode == 2) {
    evtgen = CreateBlastWaveEventGenerator(run_parameters);
  }
  evtgen->CheckSetParameters();
  cout << "Initialization complete!" << endl;

  ThermalParticleSystem* TPS = evtgen->ThermalModel()->TPS();

  // Measure time
  double wt1 = get_wall_time();

  bool infinite_mode = (run_parameters.nevents < 0);

  cout << endl;
  if (infinite_mode)
    cout << "Sampling events indefinitely (nevents < 0)..." << endl;
  else
    cout << "Sampling " << run_parameters.nevents << " events..." << endl;
  cout << "pT cuts: " << pTmin << " < pT < " << pTmax << " GeV/c" << endl;

  // Prepare statistics
  EventsProcessorSAM3 nstats(nsubs, dY, pTmin, pTmax);

  int decays = lround(run_parameters.parameters["decays"]);

  // Event loop (runs indefinitely if nevents < 0)
  for (long long event_number = 0; infinite_mode || event_number < run_parameters.nevents; ++event_number) {
    // Get primordial event (no decays yet)
    SimpleEvent evt = evtgen->GetEvent(false);

    // ProcessEvent handles decays internally:
    // 4pi charges computed before decays, accepted p/pbar after decays
    nstats.ProcessEvent(evt, *TPS, decays);

    if (infinite_mode) {
      if ((event_number + 1) % 1000 == 0) {
        cout << (event_number + 1) << " ";
        cout.flush();

        WriteToFile(prefix, nstats);
        WriteSAM3CorrectedFile(prefix, nstats);
      }
    }
    else if (run_parameters.nevents < 100
        || (event_number + 1) % (run_parameters.nevents / 100) == 0
        || (event_number + 1) % 1000 == 0) {
      cout << (event_number + 1) << " ";
      cout.flush();

      WriteToFile(prefix, nstats);
      WriteSAM3CorrectedFile(prefix, nstats);
    }
  }
  cout << endl;

  // Final write
  WriteToFile(prefix, nstats);
  WriteSAM3CorrectedFile(prefix, nstats);

  // Cleanup
  delete evtgen;
  delete TPS;

  double wt2 = get_wall_time();
  cout << "Time per single event: " << (wt2 - wt1) / run_parameters.nevents * 1.e3 << " ms" << endl;

  return 0;
}
