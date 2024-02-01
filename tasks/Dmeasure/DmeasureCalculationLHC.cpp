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

// Statistics for the D-measure


// Nch acceptances
struct NchAcceptance {
  double pTmin;
  double pTmax;
  int eta_iters;
  double dEta;
};

int eta_iters = 200;
double dEta = 0.025;

vector<NchAcceptance> nch_acceptances = {
        {0., 100., eta_iters, dEta},
        {0.2, 2.0, eta_iters, dEta},
        {0.2, 5.0, eta_iters, dEta},
        {0.6, 2.0, eta_iters, dEta},
        {0.6, 5.0, eta_iters, dEta}
};


// Process a single event
void ProcessEventForDmeasure(const SimpleEvent& event, SampleMoments::NumberStatistics& statsNchtot, vector<vector<SampleMoments::TwoNumberStatistics>>& allstatsNetSum) {
  int Nchtot = 0;
  vector<vector<int>> Nps(allstatsNetSum.size(),vector<int>(allstatsNetSum[0].size(),0));
  vector<vector<int>> Nms(allstatsNetSum.size(),vector<int>(allstatsNetSum[0].size(),0));

  // Iterate over all particles in the event
  for (const auto& particle : event.Particles) {
    // Check if the particle is a charged hadron
    int Q = parts_pdg.ParticleByPDG(particle.PDGID).ElectricCharge();
    if (Q != 0) {
      // Add the particle to the total number of charged particles
      Nchtot += 1;

      // Iterate over all acceptances
      for (int i = 0; i < nch_acceptances.size(); i++) {
        // Check if the particle is in the acceptance
        if (particle.GetPt() > nch_acceptances[i].pTmin && particle.GetPt() < nch_acceptances[i].pTmax) {
          int tindeta = abs(particle.GetEta()) / nch_acceptances[i].dEta;
          if (tindeta < Nps[i].size()) {
            // Add the particle to the number of charged particles in the acceptance
            if (Q > 0) {
              Nps[i][tindeta] += 1;
            } else {
              Nms[i][tindeta] += 1;
            }
          }
        }
      }
    }
  }
  // Compute prefix sums
  for (int i = 0; i < nch_acceptances.size(); i++) {
    for (int j = 1; j < Nps[i].size(); j++) {
      Nps[i][j] += Nps[i][j-1];
      Nms[i][j] += Nms[i][j-1];
    }
  }

  // Add observables to the statistics
  statsNchtot.AddObservation(Nchtot);
  for (int i = 0; i < nch_acceptances.size(); i++) {
    for (int j = 0; j < Nps[i].size(); j++) {
      allstatsNetSum[i][j].AddObservation(Nps[i][j] - Nms[i][j], Nps[i][j] + Nms[i][j]);
    }
  }
}

void WriteDmeasureResultsSingle(ofstream& file, SampleMoments::NumberStatistics& statsNchtot, vector<SampleMoments::TwoNumberStatistics>& statsNetSum, const NchAcceptance& nch_acceptance) {
  file << "# pT acceptance: " << nch_acceptance.pTmin << " - " << nch_acceptance.pTmax << " GeV/c" << endl;
  file << setw(15) << "dEta_acc" << " "
       << setw(15) << "D" << " "
       << setw(15) << "error" << " "
       << setw(15) << "alphach" << " "
       << setw(15) << "D'" << " "
       << setw(15) << "error" << " "
       << setw(15) << "D''" << " "
       << setw(15) << "error" << " ";
  file << endl;

  // Compute D = 4 dQ^2 / Nch, for LHC only!
  for(int ieta = 0; ieta < statsNetSum.size(); ieta++) {
    double dEtaAcc = 2. * nch_acceptance.dEta * (ieta + 1);
    file << setw(15) << dEtaAcc << " ";
    double D = 4. * statsNetSum[ieta].GetJointCumulantRatio(2, 0, 0, 1);
    double Derror = 4. * statsNetSum[ieta].GetJointCumulantRatioError(2, 0, 0, 1);
    file << setw(15) << D << " ";
    file << setw(15) << Derror << " ";
    double alphach = statsNetSum[ieta].GetMean2() / statsNchtot.GetMean();
    file << setw(15) << alphach << " ";
    double Dprime = D + 4. * alphach;
    double Dprimeerror = Derror;
    file << setw(15) << Dprime << " ";
    file << setw(15) << Dprimeerror << " ";
    double Dprimeprime = D / (1. - alphach);
    double Dprimeprimeerror = Derror / (1. - alphach);
    file << setw(15) << Dprimeprime << " ";
    file << setw(15) << Dprimeprimeerror << " ";
    file << endl;
  }
}

// Write D-measure results to file
void WriteDmeasureResults(const string& filename, SampleMoments::NumberStatistics& statsNchtot, vector<SampleMoments::TwoNumberStatistics>& statsNetSum, const NchAcceptance& nch_acceptance) {
  ofstream file(filename);
  file << "# D-measure results" << endl;
  file << "# Events: " << statsNchtot.GetNumberOfObservations() << endl;
  WriteDmeasureResultsSingle(file, statsNchtot, statsNetSum, nch_acceptance);
  file.close();
}

// Write D-measure results to file
void WriteDmeasureResultsAll(const string& filename, SampleMoments::NumberStatistics& statsNchtot, vector<vector<SampleMoments::TwoNumberStatistics>>& allstatsNetSum, const vector<NchAcceptance>& nch_acceptances) {
  ofstream file(filename);
  file << "# D-measure results" << endl;
  file << "# Events: " << statsNchtot.GetNumberOfObservations() << endl;
  for (int i = 0; i < nch_acceptances.size(); i++) {
    WriteDmeasureResultsSingle(file, statsNchtot, allstatsNetSum[i], nch_acceptances[i]);
    file << endl << endl;
  }
  file.close();
}

int main(int argc, char* argv[]) {

  cout << "Running FIST sampler version " << FistSampler_VERSION_MAJOR << "." << FistSampler_VERSION_MINOR << endl << endl;

  FistSamplerParameters run_parameters;


  string fileinput = std::string(FistSampler_INPUT_FOLDER) + "/input.ALICE.PbPb.2760.C0-5.EVHRG";
  //fileinput = std::string(FistSampler_INPUT_FOLDER) + "/input.AuAu.7.7.C0-5.EVHRG";
  //fileinput = std::string(FistSampler_INPUT_FOLDER) + "/input.HADES.AuAu.2.4.C0-5.EVHRG";
  //fileinput = std::string(FistSampler_INPUT_FOLDER) + "/input.AuAu.7.7.C70-80.EVHRG";

  if (argc > 1) {
    fileinput = string(argv[1]);
    // run_parameters.ReadParametersFromFile(fileinput);
  }

  run_parameters.ReadParametersFromFile(fileinput);
  //run_parameters.hypersurface_file = std::string(FistSampler_INPUT_FOLDER) + "/hydro/AuAu7.7/C70-80/surface_eps_0.26.dat";


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
  // Nch in full space
  SampleMoments::NumberStatistics statsNchtot;

  // Net charge in acceptance regions
  vector<vector<SampleMoments::TwoNumberStatistics>> allstatsNetSum(nch_acceptances.size());
  for(int iacc = 0; iacc < nch_acceptances.size(); iacc++) {
    allstatsNetSum[iacc].resize(nch_acceptances[iacc].eta_iters);
  }

  // Same but for primordial hadrons
  SampleMoments::NumberStatistics statsNchtotPrim;
  vector<vector<SampleMoments::TwoNumberStatistics>> allstatsNetSumPrim(nch_acceptances.size());
  for(int iacc = 0; iacc < nch_acceptances.size(); iacc++) {
    allstatsNetSumPrim[iacc].resize(nch_acceptances[iacc].eta_iters);
  }

  // Loop through the events
  for (long long event_number = 0; event_number < run_parameters.nevents; ++event_number) {
    // Sample the primordial hadrons
    SimpleEvent evt = evtgen->GetEvent(false);

    // Gather statistics for primordial hadrons
    ProcessEventForDmeasure(evt, statsNchtotPrim, allstatsNetSumPrim);

    // Perform the decays, if necessary
    if (lround(run_parameters.parameters["decays"]) != 0) {
      evt = EventGeneratorBase::PerformDecays(evt, TPS);
    }

    // Write the event to file
    if (event_writer != NULL) {
      event_writer->WriteEvent(evt);
    }

    // Process event for the D-measure
    ProcessEventForDmeasure(evt, statsNchtot, allstatsNetSum);

    // Periodically print the number of processed events on screen (every 1% or every 1000 events)
    if (run_parameters.nevents < 100 || (event_number + 1) % (run_parameters.nevents / 100) == 0 || (event_number + 1) % 1000 == 0) {
      cout << (event_number + 1) << " ";
      cout.flush();

//      for(int iacc = 0; iacc < nch_acceptances.size(); iacc++) {
//        string filename = "nchacc." + to_string(iacc) + "." + run_parameters.output_file;
//        WriteDmeasureResults(filename, statsNchtot, allstatsNetSum[iacc], nch_acceptances[iacc]);
//      }
      WriteDmeasureResultsAll("Primordial." + run_parameters.output_file, statsNchtotPrim, allstatsNetSumPrim, nch_acceptances);
      WriteDmeasureResultsAll(run_parameters.output_file, statsNchtot, allstatsNetSum, nch_acceptances);
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
