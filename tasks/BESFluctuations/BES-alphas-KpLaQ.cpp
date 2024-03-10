#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdio>
#include <set>
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
//using namespace sample_moments;

ThermalParticleSystem parts_pdg(ThermalFIST_DEFAULT_LIST_FILE);

struct AllStats {
  // Protons
  double pTminprot = 0.4, pTmaxprot = 2.0;
  double ycutprot = 0.5;
  SampleMoments::TwoNumberStatistics statspB;
  SampleMoments::TwoNumberStatistics statsppB;
  SampleMoments::TwoNumberStatistics statspmB;

  // Lambdas
  double pTminlam = 0.9, pTmaxlam = 2.0;
  double ycutL = 0.5;
  SampleMoments::TwoNumberStatistics statsLB;
  SampleMoments::TwoNumberStatistics statsLpB;
  SampleMoments::TwoNumberStatistics statsLmB;
  SampleMoments::TwoNumberStatistics statsLS;
  SampleMoments::TwoNumberStatistics statsLpS;
  SampleMoments::TwoNumberStatistics statsLmS;
  SampleMoments::TwoNumberStatistics statsLNS;
  SampleMoments::TwoNumberStatistics statsLpNS;
  SampleMoments::TwoNumberStatistics statsLmNS;

  // Kaons
  double pTminka = 0.2, pTmaxka = 1.6;
  double ycutka = 0.5;
  SampleMoments::TwoNumberStatistics statsKS;
  SampleMoments::TwoNumberStatistics statsKpS;
  SampleMoments::TwoNumberStatistics statsKmS;
  SampleMoments::TwoNumberStatistics statsKNS;
  SampleMoments::TwoNumberStatistics statsKpNS;
  SampleMoments::TwoNumberStatistics statsKmNS;

  // Net-charge
  double pTch = 0.2, pTmaxch = 2.0;
  double etacutch = 0.5;
  SampleMoments::TwoNumberStatistics statsNchQ;
  SampleMoments::TwoNumberStatistics statsNchpQ;
  SampleMoments::TwoNumberStatistics statsNchmQ;
  SampleMoments::TwoNumberStatistics statsNchNch;
  SampleMoments::TwoNumberStatistics statsNchpNch;
  SampleMoments::TwoNumberStatistics statsNchmNch;

  void ProcessEvent(const SimpleEvent& event, ThermalParticleSystem *TPS) {
    // Global charges
    int Bp = 0, Bm = 0;
    int NSp = 0, NSm = 0;
    int Sp = 0, Sm = 0;
    int NQp = 0, NQm = 0;
    int Qp = 0, Qm = 0;

    for (const auto& particle : event.Particles) {
      int B = TPS->ParticleByPDG(particle.PDGID).BaryonCharge();
      int Q = TPS->ParticleByPDG(particle.PDGID).ElectricCharge();
      int S = TPS->ParticleByPDG(particle.PDGID).Strangeness();
      if (B > 0)
        Bp += B;
      if (B < 0)
        Bm += abs(B);

      if (Q > 0) {
        Qp += Q;
        NQp += 1;
      }
      if (Q < 0) {
        Qm += abs(Q);
        NQm += 1;
      }

      if (S > 0) {
        Sp += S;
        NSp += 1;
      }
      if (S < 0) {
        Sm += abs(S);
        NSm += 1;
      }
    }

    // Protons
    int Nprot = 0, Naprot = 0;
    for (const auto& particle : event.Particles) {
      if (particle.PDGID == 2212) {
        if (particle.GetPt() > pTminprot && particle.GetPt() < pTmaxprot && abs(particle.GetY()) < ycutprot) {
          Nprot += 1;
        }
      }
      if (particle.PDGID == -2212) {
        if (particle.GetPt() > pTminprot && particle.GetPt() < pTmaxprot && abs(particle.GetY()) < ycutprot) {
          Naprot += 1;
        }
      }
    }

    // Lambda's
    int NL = 0, NAL = 0;
    for (const auto& particle : event.AllParticles) {
      if (particle.PDGID == 3122) {
        if (particle.GetPt() > pTminlam && particle.GetPt() < pTmaxlam && abs(particle.GetY()) < ycutL) {
          NL += 1;
        }
      }
      if (particle.PDGID == -3122) {
        if (particle.GetPt() > pTminlam && particle.GetPt() < pTmaxlam && abs(particle.GetY()) < ycutL) {
          NAL += 1;
        }
      }
    }

    // Kaons
    int NK = 0, NAK = 0;
    for (const auto& particle : event.AllParticles) {
      if (particle.PDGID == 321) {
        if (particle.GetPt() > pTminka && particle.GetPt() < pTmaxka && abs(particle.GetY()) < ycutka) {
          NK += 1;
        }
      }
      if (particle.PDGID == -321) {
        if (particle.GetPt() > pTminka && particle.GetPt() < pTmaxka && abs(particle.GetY()) < ycutka) {
          NAK += 1;
        }
      }
    }

    // Charged particles
    int Nchp = 0, Nchm = 0;
    int Qpacc = 0, Qmacc = 0;
    for (const auto& particle : event.Particles) {
      int Q = TPS->ParticleByPDG(particle.PDGID).ElectricCharge();
      if (Q != 0) {
        if (particle.GetPt() > pTch && particle.GetPt() < pTmaxch && abs(particle.GetEta()) < etacutch) {
          if (Q > 0) {
            Nchp += 1;
            Qpacc += Q;
          } else {
            Nchm += 1;
            Qmacc += abs(Q);
          }
        }
      }
    }

    // Fill statistics
    statspB.AddObservation(Nprot + Naprot, Bp + Bm);
    statsppB.AddObservation(Nprot, Bp + Bm);
    statspmB.AddObservation(Naprot, Bp + Bm);

    statsLB.AddObservation(NL + NAL, Bp + Bm);
    statsLpB.AddObservation(NL, Bp + Bm);
    statsLmB.AddObservation(NAL, Bp + Bm);
    statsLS.AddObservation(NL + NAL, Sp + Sm);
    statsLpS.AddObservation(NL, Sp + Sm);
    statsLmS.AddObservation(NAL, Sp + Sm);
    statsLNS.AddObservation(NL + NAL, NSp + NSm);
    statsLpNS.AddObservation(NL, NSp + NSm);
    statsLmNS.AddObservation(NAL, NSp + NSm);

    statsKS.AddObservation(NK + NAK, Sp + Sm);
    statsKpS.AddObservation(NK, Sp + Sm);
    statsKmS.AddObservation(NAK, Sp + Sm);
    statsKNS.AddObservation(NK + NAK, NSp + NSm);
    statsKpNS.AddObservation(NK, NSp + NSm);
    statsKmNS.AddObservation(NAK, NSp + NSm);

    statsNchQ.AddObservation(Qpacc + Qmacc, Qp + Qm);
    statsNchpQ.AddObservation(Qpacc, Qp + Qm);
    statsNchmQ.AddObservation(Qmacc, Qp + Qm);
    statsNchNch.AddObservation(Qpacc + Qmacc, NQp + NQm);
    statsNchpNch.AddObservation(Qpacc, NQp + NQm);
    statsNchmNch.AddObservation(Qmacc, NQp + NQm);

  }

};

FistSamplerParameters run_parameters;

void WriteAlphaResults(ostream *out, AllStats& stats) {
  ostream& file = *out;
  file << "# protons" << endl;
  file << "# pT acceptance: " << stats.pTminprot << " < pT < " << stats.pTmaxprot << " GeV/c" << endl;
  file << "# y acceptance: |y| < " << stats.ycutprot << endl;
  file << setw(15) << "ecm[GeV]" << " "
       << setw(15) << "Nprot" << " "
       << setw(15) << "error" << " "
       << setw(15) << "Naprot" << " "
       << setw(15) << "error" << " "
//       << setw(15) << "Bp" << " "
//       << setw(15) << "error" << " "
//       << setw(15) << "Bm" << " "
//       << setw(15) << "error" << " "
       << setw(15) << "alphanetp" << " "
       << setw(15) << "error" << " "
       << setw(15) << "alphap" << " "
       << setw(15) << "error" << " "
       << setw(15) << "alphaap" << " "
       << setw(15) << "error" << " ";
  file << endl;

  file << setw(15) << run_parameters.parameters["ecm"] << " ";
  file << setw(15) << stats.statsppB.GetMean1() << " ";
  file << setw(15) << stats.statsppB.GetMean1Error() << " ";
  file << setw(15) << stats.statspmB.GetMean1() << " ";
  file << setw(15) << stats.statspmB.GetMean1Error() << " ";
  file << setw(15) << stats.statspB.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statspB.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsppB.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsppB.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statspmB.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statspmB.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << endl;
  file << endl;

  file << "# lambdas" << endl;
  file << "# pT acceptance: " << stats.pTminlam << " < pT < " << stats.pTmaxlam << " GeV/c" << endl;
  file << "# y acceptance: |y| < " << stats.ycutL << endl;
  file << setw(15) << "ecm[GeV]" << " "
       << setw(15) << "NL" << " "
       << setw(15) << "error" << " "
       << setw(15) << "NAL" << " "
       << setw(15) << "error" << " ";
  file << setw(15) << "alphaBnetL" << " ";
  file << setw(15) << "error" << " ";
  file << setw(15) << "alphaBL" << " ";
  file << setw(15) << "error" << " ";
  file << setw(15) << "alphaBaL" << " ";
  file << setw(15) << "error" << " ";
  file << setw(15) << "alphaSnetL" << " ";
  file << setw(15) << "error" << " ";
  file << setw(15) << "alphaSL" << " ";
  file << setw(15) << "error" << " ";
  file << setw(15) << "alphaSaL" << " ";
  file << setw(15) << "error" << " ";
  file << endl;

  file << setw(15) << run_parameters.parameters["ecm"] << " ";
  file << setw(15) << stats.statsLpB.GetMean1() << " ";
  file << setw(15) << stats.statsLpB.GetMean1Error() << " ";
  file << setw(15) << stats.statsLmB.GetMean1() << " ";
  file << setw(15) << stats.statsLmB.GetMean1Error() << " ";
  file << setw(15) << stats.statsLB.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsLB.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsLpB.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsLpB.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsLmB.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsLmB.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsLS.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsLS.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsLpS.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsLpS.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsLmS.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsLmS.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << endl;
  file << endl;

  file << "# kaons" << endl;
  file << "# pT acceptance: " << stats.pTminka << " < pT < " << stats.pTmaxka << " GeV/c" << endl;
  file << "# y acceptance: |y| < " << stats.ycutka << endl;
  file << setw(15) << "ecm[GeV]" << " "
       << setw(15) << "NK" << " "
       << setw(15) << "error" << " "
       << setw(15) << "NAK" << " "
       << setw(15) << "error" << " ";
  file << setw(15) << "alphaSnetK" << " ";
  file << setw(15) << "error" << " ";
  file << setw(15) << "alphaSK" << " ";
  file << setw(15) << "error" << " ";
  file << setw(15) << "alphaSaK" << " ";
  file << setw(15) << "error" << " ";
  file << endl;

  file << setw(15) << run_parameters.parameters["ecm"] << " ";
  file << setw(15) << stats.statsKpS.GetMean1() << " ";
  file << setw(15) << stats.statsKpS.GetMean1Error() << " ";
  file << setw(15) << stats.statsKmS.GetMean1() << " ";
  file << setw(15) << stats.statsKmS.GetMean1Error() << " ";
  file << setw(15) << stats.statsKS.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsKS.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsKpS.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsKpS.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsKmS.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsKmS.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << endl;
  file << endl;

  file << "# net-charge" << endl;
  file << "# pT acceptance: " << stats.pTch << " < pT < " << stats.pTmaxch << " GeV/c" << endl;
  file << "# eta acceptance: |eta| < " << stats.etacutch << endl;
  file << setw(15) << "ecm[GeV]" << " "
       << setw(15) << "Qpacc" << " "
       << setw(15) << "error" << " "
       << setw(15) << "Qmacc" << " "
       << setw(15) << "error" << " ";
  file << setw(15) << "alphaQQ" << " ";
  file << setw(15) << "error" << " ";
  file << setw(15) << "alphaQpQ" << " ";
  file << setw(15) << "error" << " ";
  file << setw(15) << "alphaQmQ" << " ";
  file << setw(15) << "error" << " ";
  file << endl;

  file << setw(15) << run_parameters.parameters["ecm"] << " ";
  file << setw(15) << stats.statsNchpQ.GetMean1() << " ";
  file << setw(15) << stats.statsNchpQ.GetMean1Error() << " ";
  file << setw(15) << stats.statsNchmQ.GetMean1() << " ";
  file << setw(15) << stats.statsNchmQ.GetMean1Error() << " ";
  file << setw(15) << stats.statsNchQ.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsNchQ.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsNchpQ.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsNchpQ.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsNchmQ.GetJointMomentRatio({1,0}, {0,1}) << " ";
  file << setw(15) << stats.statsNchmQ.GetJointMomentRatioError({1,0}, {0,1}) << " ";
  file << endl;

  file << endl;
}

// Write results to file
void WriteAlphaResults(const string& filename, AllStats& stats, int nevents) {
  ofstream file(filename);
  file << "# Alpha parameter results" << endl;
  file << "# Events: " << nevents << endl;
  WriteAlphaResults(&file, stats);
  file.close();
}


int main(int argc, char* argv[]) {

  cout << "Running FIST sampler version " << FistSampler_VERSION_MAJOR << "." << FistSampler_VERSION_MINOR << endl << endl;

  string fileinput = std::string(FistSampler_INPUT_FOLDER) + "/input.AuAu.7.7.C70-80";

  run_parameters.parameters["ecm"] = 7.7;

  if (argc > 1) {
    fileinput = string(argv[1]);
    // run_parameters.ReadParametersFromFile(fileinput);
  }

  run_parameters.ReadParametersFromFile(fileinput);
  //run_parameters.hypersurface_file = std::string(FistSampler_INPUT_FOLDER) + "/hydro/AuAu7.7/C70-80/surface_eps_0.26.dat";

  // Force strong + EM decays
  run_parameters.parameters["decays"] = 3;

  if (argc > 2) {
    run_parameters.output_file = string(argv[2]);
  }

  // Output the values of all the parameters used
  run_parameters.OutputParameters();

  // Set the random seed
  RandomGenerators::SetSeed(run_parameters.randomseed);

  int fist_sampler_mode = lround(run_parameters.parameters["fist_sampler_mode"]);
  if (fist_sampler_mode < 0 || fist_sampler_mode > 2) {
    std::cout << "fist_sampler_mode of " << fist_sampler_mode << " is unsupported! " << "Aborting..." << "\n";
    exit(1);
  }

  // Cooper-Frye hypersurface. Not used for fist_sampler_mode == 2 (blast-wave)
  ParticlizationHypersurface hypersurface;

  if (fist_sampler_mode == 0) {
    // Read the Cooper-Frye hypersurface from file
    ReadHypersurfaceFromFile(run_parameters, hypersurface);

    // Check if the hypersurface is non-empty
    if (hypersurface.size() == 0) {
      std::cout << "Empty hypersurface! Aborting..." << "\n";
      exit(1);
    }
  }
  else if (fist_sampler_mode == 1) {
    // Create a Siemens-Rasmussen-Hubble hypersurface based on parameters provided in run_parameters
    CreateSiemensRasmussenHubbleHypersurface(run_parameters, hypersurface);
  }


  cout << "Initializing event generator..." << "\n";
  EventGeneratorBase* evtgen;
  
  if (fist_sampler_mode == 0 || fist_sampler_mode == 1) {
    evtgen = CreateEventGeneratorFromHypersurface(run_parameters, hypersurface);
  } 
  else if (fist_sampler_mode == 2) {
    evtgen = CreateBlastWaveEventGenerator(run_parameters);
  }
  evtgen->CheckSetParameters();
  cout << "Initialization complete!" << "\n";

  // Pointer to the particle list
  ThermalParticleSystem* TPS = evtgen->ThermalModel()->TPS();

  // Prepare the event output to file
  EventWriter* event_writer = NULL;
  if (lround(run_parameters.parameters["output_format"]) == 0)
    event_writer = new EventWriter(run_parameters.output_file);
  else if (lround(run_parameters.parameters["output_format"]) == 1)
    event_writer = new EventWriterForUrqmd(run_parameters.output_file);

  ofstream fout_events;
  if (event_writer != NULL) {
    fout_events.open(run_parameters.output_file);
  }

  // Measure time
  double wt1 = get_wall_time();

  cout << "\n";
  cout << "Sampling " << run_parameters.nevents << " events..." << endl;

  // Prepare statistics
  AllStats stats;

  // Loop through the events
  for (long long event_number = 0; event_number < run_parameters.nevents; ++event_number) {
    // Sample the primordial hadrons
    SimpleEvent evt = evtgen->GetEvent(false);

    // Perform the decays, if necessary
    if (lround(run_parameters.parameters["decays"]) != 0) {
      evt = EventGeneratorBase::PerformDecays(evt, TPS);
    }

    stats.ProcessEvent(evt, TPS);

    // Periodically print the number of processed events on screen (every 1% or every 1000 events)
    if (run_parameters.nevents < 100 || (event_number + 1) % (run_parameters.nevents / 100) == 0 || (event_number + 1) % 1000 == 0) {
      cout << (event_number + 1) << " ";
      cout.flush();

      WriteAlphaResults(&std::cout, stats);
      WriteAlphaResults(run_parameters.output_file, stats, event_number + 1);
    }
  }
  cout << endl;

  // Cleanup
  delete evtgen;
  delete TPS;

  // Time performace
  double wt2 = get_wall_time();

  cout << "Time per single event: " << (wt2 - wt1) / run_parameters.nevents * 1.e3 << " ms" << "\n";

  return 0;
}
