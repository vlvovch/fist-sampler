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

  // Joint statistics of net-proton X = N_p - N_pbar and 4pi baryon number B,
  // used both for direct-MC net-proton cumulants and for the SAM-3.0
  // B-canonical correction to κ_n[X].
  vector<TwoNumberStatistics> statsXB;

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
    statsXB.resize(nbins);
  }


  /// Per-event summary: prefix-summed Np/Npbar over rapidity bins and 4pi charges.
  /// Shared across multiple processors (e.g. jackknife blocks) so that the
  /// expensive decay step is done only once per event.
  struct EventData {
    vector<int> Np;
    vector<int> Npbar;
    int Btot = 0;
    int Qtot = 0;
    int Stot = 0;
  };

  /// Expensive step: perform decays, bin accepted p/pbar, sum 4pi charges
  /// (computed before decays so weak decays don't leak S conservation).
  static EventData ComputeEventData(const SimpleEvent& evt, ThermalParticleSystem& TPS,
                                    int decays, int nsubs, double dY,
                                    double pTmin, double pTmax) {
    int nbins = nsubs / 2;
    EventData ed;
    ed.Np.assign(nbins, 0);
    ed.Npbar.assign(nbins, 0);

    for (const SimpleParticle& part : evt.Particles) {
      const ThermalParticle& prop = TPS.ParticleByPDG(part.PDGID);
      ed.Btot += prop.BaryonCharge();
      ed.Qtot += prop.ElectricCharge();
      ed.Stot += prop.Strangeness();
    }

    SimpleEvent evt_final = evt;
    if (decays != 0) {
      evt_final = EventGeneratorBase::PerformDecays(evt, &TPS);
    }

    for (const SimpleParticle& part : evt_final.Particles) {
      double pT = part.GetPt();
      if (pT >= pTmin && pT <= pTmax) {
        int tindY = floor(abs(part.GetY()) / dY);
        if (tindY >= 0 && tindY < nbins) {
          if (part.PDGID == 2212)
            ed.Np[tindY]++;
          if (part.PDGID == -2212)
            ed.Npbar[tindY]++;
        }
      }
    }

    for (int i = 1; i < nbins; ++i) {
      ed.Np[i]    += ed.Np[i - 1];
      ed.Npbar[i] += ed.Npbar[i - 1];
    }
    return ed;
  }

  /// Cheap step: fold a pre-computed event into this processor's statistics.
  void AddEventData(const EventData& ed) {
    nevents++;
    int nbins = m_nsubs / 2;
    for (int isub = 0; isub < nbins; ++isub) {
      statsPPbar[isub].AddObservation(ed.Np[isub], ed.Npbar[isub]);
      statsPB[isub].AddObservation(ed.Np[isub], ed.Btot);
      statsPQ[isub].AddObservation(ed.Np[isub], ed.Qtot);
      statsPS[isub].AddObservation(ed.Np[isub], ed.Stot);
      statsPbarB[isub].AddObservation(ed.Npbar[isub], ed.Btot);
      statsPbarQ[isub].AddObservation(ed.Npbar[isub], ed.Qtot);
      statsPbarS[isub].AddObservation(ed.Npbar[isub], ed.Stot);
      statsXB   [isub].AddObservation(ed.Np[isub] - ed.Npbar[isub], ed.Btot);
    }
    statsBQ.AddObservation(ed.Btot, ed.Qtot);
    statsBS.AddObservation(ed.Btot, ed.Stot);
    statsQS.AddObservation(ed.Qtot, ed.Stot);
    statsB.AddObservation(ed.Btot);
    statsQ.AddObservation(ed.Qtot);
    statsS.AddObservation(ed.Stot);
  }

  /// Backward-compatible single-call interface.
  void ProcessEvent(const SimpleEvent& evt, ThermalParticleSystem& TPS, int decays) {
    EventData ed = ComputeEventData(evt, TPS, decays, m_nsubs, m_dY, m_pTmin, m_pTmax);
    AddEventData(ed);
  }
};


// ============================================================
// SAM-3.0 value container — one entry per ensemble per y-bin.
// Ensembles (index 0..6): gce, Bcan, Qcan, Scan, BQcan, BScan, BQScan.
// Used both for central values and for per-block jackknife samples.
// ============================================================
constexpr int N_SAM3_ENSEMBLES = 7;
inline const string& SAM3EnsembleName(int i) {
  static const string names[N_SAM3_ENSEMBLES] = {
    "gce", "Bcan", "Qcan", "Scan", "BQcan", "BScan", "BQScan"
  };
  return names[i];
}

struct SAM3Cumulants {
  int nbins = 0;
  vector<double> meanNp, meanNpb;               // [nbins]
  vector<vector<double>> k2p, k2pb, k11;        // [7][nbins]
  // Direct-MC cumulant ratios κ_3/κ_1 and κ_4/κ_2 for p and p̄.
  // Not SAM-3.0 corrected; per-bin only.
  vector<double> r_k3k1_p, r_k4k2_p;            // [nbins]
  vector<double> r_k3k1_pb, r_k4k2_pb;          // [nbins]
  // Direct-MC net-proton ratios for X = N_p − N_pbar.
  vector<double> r_k2X_skell, r_k3X_k1X, r_k4X_k2X;  // [nbins]
};

// Compute all SAM-3.0 corrected 2nd-order cumulants for a given processor.
// Returns an empty (nbins=0) result if fewer than 2 events have been added.
static SAM3Cumulants ComputeSAM3Cumulants(EventsProcessorSAM3& stats) {
  SAM3Cumulants out;
  if (stats.nevents < 2) return out;

  int nbins = stats.m_nsubs / 2;
  out.nbins = nbins;
  out.meanNp.assign(nbins, 0.);
  out.meanNpb.assign(nbins, 0.);
  out.k2p .assign(N_SAM3_ENSEMBLES, vector<double>(nbins, 0.));
  out.k2pb.assign(N_SAM3_ENSEMBLES, vector<double>(nbins, 0.));
  out.k11 .assign(N_SAM3_ENSEMBLES, vector<double>(nbins, 0.));
  out.r_k3k1_p .assign(nbins, 0.);
  out.r_k4k2_p .assign(nbins, 0.);
  out.r_k3k1_pb.assign(nbins, 0.);
  out.r_k4k2_pb.assign(nbins, 0.);
  out.r_k2X_skell.assign(nbins, 0.);
  out.r_k3X_k1X  .assign(nbins, 0.);
  out.r_k4X_k2X  .assign(nbins, 0.);

  double varB  = stats.statsB.GetCentralMoment(2);
  double varQ  = stats.statsQ.GetCentralMoment(2);
  double varS  = stats.statsS.GetCentralMoment(2);
  double covBQ = stats.statsBQ.GetJointCentralMoment(1, 1);
  double covBS = stats.statsBS.GetJointCentralMoment(1, 1);
  double covQS = stats.statsQS.GetJointCentralMoment(1, 1);

  double K[3][3] = {
    {varB,  covBQ, covBS},
    {covBQ, varQ,  covQS},
    {covBS, covQS, varS}
  };

  auto invertSubmatrix = [&K](const vector<int>& idx) -> vector<double> {
    int N = (int)idx.size();
    if (N == 1) return {1.0 / K[idx[0]][idx[0]]};
    if (N == 2) {
      double a = K[idx[0]][idx[0]], b = K[idx[0]][idx[1]];
      double d = K[idx[1]][idx[1]];
      double det = a * d - b * b;
      return {d / det, -b / det, -b / det, a / det};
    }
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
  };

  // Scenarios in the order matching SAM3EnsembleName() (indices 1..6).
  vector<vector<int>> scen_idx = {
    {0}, {1}, {2},
    {0, 1}, {0, 2},
    {0, 1, 2}
  };
  vector<vector<double>> Kinvs(scen_idx.size());
  for (size_t i = 0; i < scen_idx.size(); ++i)
    Kinvs[i] = invertSubmatrix(scen_idx[i]);

  auto correction = [](const vector<double>& Kinv, const vector<int>& idx,
                        const double cpi[3], const double cpj[3]) -> double {
    int N = (int)idx.size();
    double r = 0.;
    for (int a = 0; a < N; ++a)
      for (int b = 0; b < N; ++b)
        r += cpi[idx[a]] * Kinv[a * N + b] * cpj[idx[b]];
    return r;
  };

  for (int isub = 0; isub < nbins; ++isub) {
    auto& ppbar = stats.statsPPbar[isub];
    out.meanNp[isub]  = ppbar.GetMean1();
    out.meanNpb[isub] = ppbar.GetMean2();

    double varNp    = ppbar.GetJointCumulant(2, 0);
    double varNpb   = ppbar.GetJointCumulant(0, 2);
    double covNpNpb = ppbar.GetJointCumulant(1, 1);

    double cpNp[3]  = {
      stats.statsPB[isub].GetJointCumulant(1, 1),
      stats.statsPQ[isub].GetJointCumulant(1, 1),
      stats.statsPS[isub].GetJointCumulant(1, 1)
    };
    double cpNpb[3] = {
      stats.statsPbarB[isub].GetJointCumulant(1, 1),
      stats.statsPbarQ[isub].GetJointCumulant(1, 1),
      stats.statsPbarS[isub].GetJointCumulant(1, 1)
    };

    // GCE (ensemble 0): no correction
    out.k2p [0][isub] = varNp;
    out.k2pb[0][isub] = varNpb;
    out.k11 [0][isub] = covNpNpb;

    // 6 SAM-3.0-corrected scenarios (ensembles 1..6)
    for (size_t iS = 0; iS < scen_idx.size(); ++iS) {
      double dVarNp    = correction(Kinvs[iS], scen_idx[iS], cpNp,  cpNp);
      double dVarNpb   = correction(Kinvs[iS], scen_idx[iS], cpNpb, cpNpb);
      double dCovNpNpb = correction(Kinvs[iS], scen_idx[iS], cpNp,  cpNpb);

      out.k2p [iS + 1][isub] = varNp    - dVarNp;
      out.k2pb[iS + 1][isub] = varNpb   - dVarNpb;
      out.k11 [iS + 1][isub] = covNpNpb - dCovNpNpb;
    }

    // Direct-MC k3/k1 and k4/k2 ratios (ensemble-independent).
    out.r_k3k1_p [isub] = ppbar.GetJointCumulantRatio(3, 0, 1, 0);
    out.r_k4k2_p [isub] = ppbar.GetJointCumulantRatio(4, 0, 2, 0);
    out.r_k3k1_pb[isub] = ppbar.GetJointCumulantRatio(0, 3, 0, 1);
    out.r_k4k2_pb[isub] = ppbar.GetJointCumulantRatio(0, 4, 0, 2);

    // Direct-MC net-proton ratios (from statsXB).
    auto& sXB = stats.statsXB[isub];
    double skellam = ppbar.GetMean1() + ppbar.GetMean2();
    out.r_k2X_skell[isub] = (skellam > 0.) ? sXB.GetJointCumulant(2, 0) / skellam : 0.;
    out.r_k3X_k1X  [isub] = sXB.GetJointCumulantRatio(3, 0, 1, 0);
    out.r_k4X_k2X  [isub] = sXB.GetJointCumulantRatio(4, 0, 2, 0);
  }
  return out;
}


// ============================================================
// SAM-3.0 B-canonical correction for 2nd-4th cumulants of a single
// observable N, derived via saddle-point from the joint (N, B) CGF at
// fixed B = <B>.  Single-charge specialization; requires κ^gce_{02}[B] ≠ 0.
//
//   A  = κ_{11}/κ_{02},   M = κ_{21} − 2 κ_{12} A + κ_{03} A²
//   κ̃_2 = κ_{20} − κ_{11}²/κ_{02}
//   κ̃_3 = κ_{30} − 3 κ_{21} A + 3 κ_{12} A² − κ_{03} A³
//   κ̃_4 = κ_{40} − 4 κ_{31} A + 6 κ_{22} A² − 4 κ_{13} A³ + κ_{04} A⁴
//          − 3 M²/κ_{02}
//
// Indexing κ_{ij}: i derivatives w.r.t. N, j w.r.t. B.
// ============================================================
struct SAM3HigherOrderBcan { double k2, k3, k4; };

static SAM3HigherOrderBcan SAM3BcanK2K3K4(
    double kN2, double kN3, double kN4,
    double kB2, double kB3, double kB4,
    double k11, double k21, double k12,
    double k31, double k13, double k22)
{
  double A  = k11 / kB2;
  double A2 = A * A, A3 = A2 * A, A4 = A2 * A2;
  double M  = k21 - 2.0 * k12 * A + kB3 * A2;
  SAM3HigherOrderBcan r;
  r.k2 = kN2 - k11 * k11 / kB2;
  r.k3 = kN3 - 3.0 * k21 * A + 3.0 * k12 * A2 - kB3 * A3;
  r.k4 = kN4 - 4.0 * k31 * A + 6.0 * k22 * A2 - 4.0 * k13 * A3 + kB4 * A4
             - 3.0 * M * M / kB2;
  return r;
}


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
       << setw(w) << "k3(Np)"       << setw(w) << "k3(Np)_err"
       << setw(w) << "k4(Np)"       << setw(w) << "k4(Np)_err"
       << setw(w) << "k3(Npb)"      << setw(w) << "k3(Npb)_err"
       << setw(w) << "k4(Npb)"      << setw(w) << "k4(Npb)_err"
       << setw(w) << "k3Np/k1Np"    << setw(w) << "k3Np/k1Np_er"
       << setw(w) << "k4Np/k2Np"    << setw(w) << "k4Np/k2Np_er"
       << setw(w) << "k3Npb/k1Npb"  << setw(w) << "k3NpbK1Npb_e"
       << setw(w) << "k4Npb/k2Npb"  << setw(w) << "k4NpbK2Npb_e"
       << setw(w) << "k2(X)"        << setw(w) << "k2(X)_err"
       << setw(w) << "k3(X)"        << setw(w) << "k3(X)_err"
       << setw(w) << "k4(X)"        << setw(w) << "k4(X)_err"
       << setw(w) << "k2X/Skellam"  << setw(w) << "k2XSkell_err"
       << setw(w) << "k3X/k1X"      << setw(w) << "k3X/k1X_err"
       << setw(w) << "k4X/k2X"      << setw(w) << "k4X/k2X_err"
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

    // ---- 3rd and 4th cumulants of N_p, N_pbar + ratios k3/k1, k4/k2 ----
    double k3Np       = ppbar.GetJointCumulant(3, 0);
    double k3NpErr    = ppbar.GetJointCumulantError(3, 0);
    double k4Np       = ppbar.GetJointCumulant(4, 0);
    double k4NpErr    = ppbar.GetJointCumulantError(4, 0);
    double k3Npb      = ppbar.GetJointCumulant(0, 3);
    double k3NpbErr   = ppbar.GetJointCumulantError(0, 3);
    double k4Npb      = ppbar.GetJointCumulant(0, 4);
    double k4NpbErr   = ppbar.GetJointCumulantError(0, 4);

    double k3k1_p     = ppbar.GetJointCumulantRatio     (3, 0, 1, 0);
    double k3k1_pErr  = ppbar.GetJointCumulantRatioError(3, 0, 1, 0);
    double k4k2_p     = ppbar.GetJointCumulantRatio     (4, 0, 2, 0);
    double k4k2_pErr  = ppbar.GetJointCumulantRatioError(4, 0, 2, 0);
    double k3k1_pb    = ppbar.GetJointCumulantRatio     (0, 3, 0, 1);
    double k3k1_pbErr = ppbar.GetJointCumulantRatioError(0, 3, 0, 1);
    double k4k2_pb    = ppbar.GetJointCumulantRatio     (0, 4, 0, 2);
    double k4k2_pbErr = ppbar.GetJointCumulantRatioError(0, 4, 0, 2);

    // Net-proton X = N_p − N_pbar: raw cumulants (from statsXB) plus three
    // ratios.  k3X/k1X and k4X/k2X use the correlation-aware ratio error;
    // Skellam = <N_p> + <N_pbar> is effectively exact (1st moments), so the
    // naive propagated error is fine.
    auto& sXB = stats.statsXB[isub];
    double k2X        = sXB.GetJointCumulant(2, 0);
    double k2XErr     = sXB.GetJointCumulantError(2, 0);
    double k3X        = sXB.GetJointCumulant(3, 0);
    double k3XErr     = sXB.GetJointCumulantError(3, 0);
    double k4X        = sXB.GetJointCumulant(4, 0);
    double k4XErr     = sXB.GetJointCumulantError(4, 0);

    double skellam      = meanNp + meanNpb;
    double k2X_skell    = (skellam > 0.) ? k2X     / skellam : 0.;
    double k2X_skellErr = (skellam > 0.) ? k2XErr  / skellam : 0.;

    double k3X_k1X    = sXB.GetJointCumulantRatio     (3, 0, 1, 0);
    double k3X_k1XErr = sXB.GetJointCumulantRatioError(3, 0, 1, 0);
    double k4X_k2X    = sXB.GetJointCumulantRatio     (4, 0, 2, 0);
    double k4X_k2XErr = sXB.GetJointCumulantRatioError(4, 0, 2, 0);

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
         << setw(w) << k3Np       << setw(w) << k3NpErr
         << setw(w) << k4Np       << setw(w) << k4NpErr
         << setw(w) << k3Npb      << setw(w) << k3NpbErr
         << setw(w) << k4Npb      << setw(w) << k4NpbErr
         << setw(w) << k3k1_p     << setw(w) << k3k1_pErr
         << setw(w) << k4k2_p     << setw(w) << k4k2_pErr
         << setw(w) << k3k1_pb    << setw(w) << k3k1_pbErr
         << setw(w) << k4k2_pb    << setw(w) << k4k2_pbErr
         << setw(w) << k2X         << setw(w) << k2XErr
         << setw(w) << k3X         << setw(w) << k3XErr
         << setw(w) << k4X         << setw(w) << k4XErr
         << setw(w) << k2X_skell   << setw(w) << k2X_skellErr
         << setw(w) << k3X_k1X     << setw(w) << k3X_k1XErr
         << setw(w) << k4X_k2X     << setw(w) << k4X_k2XErr
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
void WriteSAM3CorrectedFile(const string& prefix, EventsProcessorSAM3& stats,
                            bool gce_mode = false) {
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
  if (gce_mode) {
    fout << setw(w) << "k3p_gce"     << setw(w) << "k4p_gce"
         << setw(w) << "k3pb_gce"    << setw(w) << "k4pb_gce"
         << setw(w) << "k3p_Bcan"    << setw(w) << "k4p_Bcan"
         << setw(w) << "k3pb_Bcan"   << setw(w) << "k4pb_Bcan"
         << setw(w) << "k2X_gce"     << setw(w) << "k3X_gce"    << setw(w) << "k4X_gce"
         << setw(w) << "k2X_Bcan"    << setw(w) << "k3X_Bcan"   << setw(w) << "k4X_Bcan";
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

    if (gce_mode) {
      // GCE cumulants of N_p, N_pbar and X = N_p − N_pbar (direct from sample).
      double k3p_gce  = ppbar.GetJointCumulant(3, 0);
      double k4p_gce  = ppbar.GetJointCumulant(4, 0);
      double k3pb_gce = ppbar.GetJointCumulant(0, 3);
      double k4pb_gce = ppbar.GetJointCumulant(0, 4);

      auto& sXB = stats.statsXB[isub];
      double k2X_gce = sXB.GetJointCumulant(2, 0);
      double k3X_gce = sXB.GetJointCumulant(3, 0);
      double k4X_gce = sXB.GetJointCumulant(4, 0);

      // Joint cumulants with 4π baryon charge B (κ_{i,j} = ∂^i_N ∂^j_B K).
      auto& sPB  = stats.statsPB [isub];
      auto& sPbB = stats.statsPbarB[isub];

      double kB2  = varB;
      double kB3  = stats.statsB.GetCumulant(3);
      double kB4  = stats.statsB.GetCumulant(4);

      double pB_11 = sPB.GetJointCumulant(1, 1);
      double pB_21 = sPB.GetJointCumulant(2, 1);
      double pB_12 = sPB.GetJointCumulant(1, 2);
      double pB_31 = sPB.GetJointCumulant(3, 1);
      double pB_13 = sPB.GetJointCumulant(1, 3);
      double pB_22 = sPB.GetJointCumulant(2, 2);

      double pbB_11 = sPbB.GetJointCumulant(1, 1);
      double pbB_21 = sPbB.GetJointCumulant(2, 1);
      double pbB_12 = sPbB.GetJointCumulant(1, 2);
      double pbB_31 = sPbB.GetJointCumulant(3, 1);
      double pbB_13 = sPbB.GetJointCumulant(1, 3);
      double pbB_22 = sPbB.GetJointCumulant(2, 2);

      double XB_11 = sXB.GetJointCumulant(1, 1);
      double XB_21 = sXB.GetJointCumulant(2, 1);
      double XB_12 = sXB.GetJointCumulant(1, 2);
      double XB_31 = sXB.GetJointCumulant(3, 1);
      double XB_13 = sXB.GetJointCumulant(1, 3);
      double XB_22 = sXB.GetJointCumulant(2, 2);

      SAM3HigherOrderBcan p_bcan = SAM3BcanK2K3K4(
        varNp, k3p_gce, k4p_gce, kB2, kB3, kB4,
        pB_11, pB_21, pB_12, pB_31, pB_13, pB_22);
      SAM3HigherOrderBcan pb_bcan = SAM3BcanK2K3K4(
        varNpb, k3pb_gce, k4pb_gce, kB2, kB3, kB4,
        pbB_11, pbB_21, pbB_12, pbB_31, pbB_13, pbB_22);
      SAM3HigherOrderBcan X_bcan = SAM3BcanK2K3K4(
        k2X_gce, k3X_gce, k4X_gce, kB2, kB3, kB4,
        XB_11, XB_21, XB_12, XB_31, XB_13, XB_22);

      fout << setw(w) << k3p_gce    << setw(w) << k4p_gce
           << setw(w) << k3pb_gce   << setw(w) << k4pb_gce
           << setw(w) << p_bcan.k3  << setw(w) << p_bcan.k4
           << setw(w) << pb_bcan.k3 << setw(w) << pb_bcan.k4
           << setw(w) << k2X_gce    << setw(w) << k3X_gce     << setw(w) << k4X_gce
           << setw(w) << X_bcan.k2  << setw(w) << X_bcan.k3   << setw(w) << X_bcan.k4;
    }
    fout << endl;
  }

  fout.close();
}


// ============================================================
// Jackknife writer: block-resample the event list to estimate the
// statistical uncertainty on each SAM-3.0 corrected cumulant, accounting
// for correlations between the input moments (they all come from the same
// events).
//
// Method: N disjoint event blocks, populated round-robin during sampling.
// Each block gives an independent SAM-3.0 estimate θ_i.  The variance of
// the full-sample estimator is estimated as
//     Var(θ_full) = Σ (θ_i − θ̄)² / (N (N − 1))
// Central values reported are from the full-sample processor (unbiased for
// non-linear estimators).
// ============================================================
void WriteSAM3JackknifeFile(const string& prefix,
                            EventsProcessorSAM3& nstats,
                            vector<EventsProcessorSAM3>& jk_blocks) {
  ofstream fout(prefix + ".SAM3-jackknife.dat");
  int w = 15;

  SAM3Cumulants central = ComputeSAM3Cumulants(nstats);
  if (central.nbins == 0) {
    fout << "# Not enough events in the full processor to compute cumulants yet." << endl;
    return;
  }
  int nbins = central.nbins;

  int N_blocks = (int)jk_blocks.size();
  vector<SAM3Cumulants> block_results(N_blocks);
  int n_good = 0;
  for (int ib = 0; ib < N_blocks; ++ib) {
    block_results[ib] = ComputeSAM3Cumulants(jk_blocks[ib]);
    if (block_results[ib].nbins == nbins) ++n_good;
  }

  fout << "# SAM-3.0 cumulants with jackknife statistical errors." << endl;
  fout << "# Error = stddev(block estimates) / sqrt(N_blocks - 1)," << endl;
  fout << "#   i.e. Var(full) = Sigma (theta_i - theta_bar)^2 / (N (N - 1))." << endl;
  fout << "# Events total: " << nstats.nevents
       << "   Blocks: " << N_blocks << " (populated: " << n_good << ")" << endl;
  fout << "# pT cuts: " << nstats.m_pTmin << " < pT < " << nstats.m_pTmax << " GeV/c" << endl;
  fout << "# Ensembles: gce, Bcan, Qcan, Scan, BQcan, BScan, BQScan" << endl;
  fout << "#" << endl;

  if (n_good < 2) {
    fout << "# Need at least 2 populated blocks; retry after more events." << endl;
    return;
  }

  // Header row.
  fout << setw(w) << "ycut"
       << setw(w) << "<Np>"
       << setw(w) << "<Npbar>";
  for (int ie = 0; ie < N_SAM3_ENSEMBLES; ++ie) {
    const string& n = SAM3EnsembleName(ie);
    fout << setw(w) << ("k2p_" + n)   << setw(w) << ("k2p_" + n + "_e")
         << setw(w) << ("k2pb_" + n)  << setw(w) << ("k2pb_" + n + "_e")
         << setw(w) << ("k11_" + n)   << setw(w) << ("k11_" + n + "_e");
  }
  fout << setw(w) << "k3p/k1p"       << setw(w) << "k3p/k1p_e"
       << setw(w) << "k4p/k2p"       << setw(w) << "k4p/k2p_e"
       << setw(w) << "k3pb/k1pb"     << setw(w) << "k3pb/k1pb_e"
       << setw(w) << "k4pb/k2pb"     << setw(w) << "k4pb/k2pb_e"
       << setw(w) << "k2X/Skellam"   << setw(w) << "k2X/Skell_e"
       << setw(w) << "k3X/k1X"       << setw(w) << "k3X/k1X_e"
       << setw(w) << "k4X/k2X"       << setw(w) << "k4X/k2X_e";
  fout << endl;

  // jk_error: std.dev of the full-sample estimate from the N block estimates.
  //   Var(full) = sample_var(block) / N_blocks
  //             = Σ(v_i - v̄)² / (N_blocks × (N_blocks - 1)).
  // Non-finite block values are skipped (e.g. ratios like k3/k1 where the
  // denominator is zero in a block with zero antiprotons at small y_cut).
  auto jk_error = [&](int ensemble, int isub, int quantity) -> double {
    double sum = 0., sum_sq = 0.;
    int n = 0;
    for (int ib = 0; ib < N_blocks; ++ib) {
      if (block_results[ib].nbins != nbins) continue;
      double v = (quantity == 0) ? block_results[ib].k2p [ensemble][isub]
               : (quantity == 1) ? block_results[ib].k2pb[ensemble][isub]
                                 : block_results[ib].k11 [ensemble][isub];
      if (!std::isfinite(v)) continue;
      sum    += v;
      sum_sq += v * v;
      ++n;
    }
    if (n < 2) return 0.;
    double mean = sum / n;
    double ss   = sum_sq - n * mean * mean;
    if (ss < 0.) ss = 0.;            // floating-point safety
    return std::sqrt(ss / ((double)n * (double)(n - 1)));
  };

  // Same formula as jk_error, applied to an arbitrary per-bin member selected
  // by `pick` (returns one double from a block's SAM3Cumulants).
  auto jk_error_ratio = [&](int isub,
                            double (*pick)(const SAM3Cumulants&, int)) -> double {
    double sum = 0., sum_sq = 0.;
    int n = 0;
    for (int ib = 0; ib < N_blocks; ++ib) {
      if (block_results[ib].nbins != nbins) continue;
      double v = pick(block_results[ib], isub);
      if (!std::isfinite(v)) continue;
      sum    += v;
      sum_sq += v * v;
      ++n;
    }
    if (n < 2) return 0.;
    double mean = sum / n;
    double ss   = sum_sq - n * mean * mean;
    if (ss < 0.) ss = 0.;
    return std::sqrt(ss / ((double)n * (double)(n - 1)));
  };

  auto pick_k3k1_p    = [](const SAM3Cumulants& c, int i) { return c.r_k3k1_p   [i]; };
  auto pick_k4k2_p    = [](const SAM3Cumulants& c, int i) { return c.r_k4k2_p   [i]; };
  auto pick_k3k1_pb   = [](const SAM3Cumulants& c, int i) { return c.r_k3k1_pb  [i]; };
  auto pick_k4k2_pb   = [](const SAM3Cumulants& c, int i) { return c.r_k4k2_pb  [i]; };
  auto pick_k2X_skell = [](const SAM3Cumulants& c, int i) { return c.r_k2X_skell[i]; };
  auto pick_k3X_k1X   = [](const SAM3Cumulants& c, int i) { return c.r_k3X_k1X  [i]; };
  auto pick_k4X_k2X   = [](const SAM3Cumulants& c, int i) { return c.r_k4X_k2X  [i]; };

  for (int isub = 0; isub < nbins; ++isub) {
    double ycut = (isub + 1) * nstats.m_dY;
    fout << setw(w) << ycut
         << setw(w) << central.meanNp[isub]
         << setw(w) << central.meanNpb[isub];
    for (int ie = 0; ie < N_SAM3_ENSEMBLES; ++ie) {
      fout << setw(w) << central.k2p [ie][isub] << setw(w) << jk_error(ie, isub, 0)
           << setw(w) << central.k2pb[ie][isub] << setw(w) << jk_error(ie, isub, 1)
           << setw(w) << central.k11 [ie][isub] << setw(w) << jk_error(ie, isub, 2);
    }
    fout << setw(w) << central.r_k3k1_p [isub]   << setw(w) << jk_error_ratio(isub, pick_k3k1_p  )
         << setw(w) << central.r_k4k2_p [isub]   << setw(w) << jk_error_ratio(isub, pick_k4k2_p  )
         << setw(w) << central.r_k3k1_pb[isub]   << setw(w) << jk_error_ratio(isub, pick_k3k1_pb )
         << setw(w) << central.r_k4k2_pb[isub]   << setw(w) << jk_error_ratio(isub, pick_k4k2_pb )
         << setw(w) << central.r_k2X_skell[isub] << setw(w) << jk_error_ratio(isub, pick_k2X_skell)
         << setw(w) << central.r_k3X_k1X  [isub] << setw(w) << jk_error_ratio(isub, pick_k3X_k1X )
         << setw(w) << central.r_k4X_k2X  [isub] << setw(w) << jk_error_ratio(isub, pick_k4X_k2X );
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

  int Bcan = lround(run_parameters.parameters["Bcanonical"]);
  int Qcan = lround(run_parameters.parameters["Qcanonical"]);
  int Scan = lround(run_parameters.parameters["Scanonical"]);
  const bool gce_mode = (!Bcan && !Qcan && !Scan);
  {
    string ensemble_suffix;
    if (gce_mode)
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

  // Prepare statistics: full-sample processor + N_JK jackknife blocks.
  EventsProcessorSAM3 nstats(nsubs, dY, pTmin, pTmax);

  const int N_JK = 20;
  vector<EventsProcessorSAM3> jk_blocks;
  jk_blocks.reserve(N_JK);
  for (int ib = 0; ib < N_JK; ++ib)
    jk_blocks.emplace_back(nsubs, dY, pTmin, pTmax);

  int decays = lround(run_parameters.parameters["decays"]);

  // Event loop (runs indefinitely if nevents < 0)
  for (long long event_number = 0; infinite_mode || event_number < run_parameters.nevents; ++event_number) {
    // Get primordial event (no decays yet).
    SimpleEvent evt = evtgen->GetEvent(false);

    // Compute the per-event summary ONCE (does the expensive decay step),
    // then fold into both the full processor and one jackknife block.
    EventsProcessorSAM3::EventData ed =
      EventsProcessorSAM3::ComputeEventData(evt, *TPS, decays, nsubs, dY, pTmin, pTmax);

    nstats.AddEventData(ed);
    jk_blocks[event_number % N_JK].AddEventData(ed);

    if (infinite_mode) {
      if ((event_number + 1) % 1000 == 0) {
        cout << (event_number + 1) << " ";
        cout.flush();

        WriteToFile(prefix, nstats);
        WriteSAM3CorrectedFile(prefix, nstats, gce_mode);
        WriteSAM3JackknifeFile(prefix, nstats, jk_blocks);
      }
    }
    else if ((event_number + 1) % 100 == 0) {
      cout << (event_number + 1) << " ";
      cout.flush();

      WriteToFile(prefix, nstats);
      WriteSAM3CorrectedFile(prefix, nstats, gce_mode);
      WriteSAM3JackknifeFile(prefix, nstats, jk_blocks);
    }
  }
  cout << endl;

  // Final write
  WriteToFile(prefix, nstats);
  WriteSAM3CorrectedFile(prefix, nstats, gce_mode);
  WriteSAM3JackknifeFile(prefix, nstats, jk_blocks);

  // Cleanup
  delete evtgen;
  delete TPS;

  double wt2 = get_wall_time();
  cout << "Time per single event: " << (wt2 - wt1) / run_parameters.nevents * 1.e3 << " ms" << endl;

  return 0;
}
