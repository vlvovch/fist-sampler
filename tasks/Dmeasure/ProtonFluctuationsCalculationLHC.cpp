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

// Nch acceptances
struct ProtonAcceptance {
  double pmin;
  double pmax;
  int eta_iters;
  double dEta;
};

int eta_iters = 200;
double dEta = 0.025;

vector<ProtonAcceptance> proton_acceptances = {
        {0.6, 1.5, eta_iters, dEta},
        {0.6, 2.0, eta_iters, dEta},
        {0.0, 1e10, eta_iters, dEta}
};


// Process a single event
void ProcessEventProtons(const SimpleEvent& event, vector<vector<SampleMoments::TwoNumberStatistics>>& allstatsNetSum) {
  vector<vector<int>> Nps(allstatsNetSum.size(),vector<int>(allstatsNetSum[0].size(),0));
  vector<vector<int>> Nms(allstatsNetSum.size(),vector<int>(allstatsNetSum[0].size(),0));

  // Iterate over all particles in the event
  for (const auto& particle : event.Particles) {
    // Check if the particle is a proton or antiproton
    if (abs(particle.PDGID) == 2212) {
      // Add the particle to the total number of charged particles
      // Iterate over all acceptances
      for (int i = 0; i < proton_acceptances.size(); i++) {
        // Check if the particle is in the acceptance
        if (particle.GetP() > proton_acceptances[i].pmin && particle.GetP() < proton_acceptances[i].pmax) {
          int tindeta = abs(particle.GetEta()) / proton_acceptances[i].dEta;
          if (tindeta < Nps[i].size()) {
            // Add the particle to the number of charged particles in the acceptance
            if (particle.PDGID > 0) {
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
  for (int i = 0; i < proton_acceptances.size(); i++) {
    for (int j = 1; j < Nps[i].size(); j++) {
      Nps[i][j] += Nps[i][j-1];
      Nms[i][j] += Nms[i][j-1];
    }
  }

  // Add observables to the statistics
  for (int i = 0; i < proton_acceptances.size(); i++) {
    for (int j = 0; j < Nps[i].size(); j++) {
      allstatsNetSum[i][j].AddObservation(Nps[i][j] - Nms[i][j], Nps[i][j] + Nms[i][j]);
    }
  }
}

void WriteProtonResultsSingle(ofstream& file, vector<SampleMoments::TwoNumberStatistics>& statsNetSum, const ProtonAcceptance& proton_acceptance) {
  file << "# p acceptance: " << proton_acceptance.pmin << " - " << proton_acceptance.pmax << " GeV/c" << endl;
  file << setw(15) << "dEta_acc" << " "
       << setw(15) << "eta_cut" << " "
       << setw(15) << "k2[net]/Sk" << " "
       << setw(15) << "error" << " "
       << setw(15) << "k2[sum]/Sk" << " "
       << setw(15) << "error" << " "
       << setw(15) << "c2net" << " "
       << setw(15) << "error" << " ";
  file << endl;

  // Compute D = 4 dQ^2 / Nch, for LHC only!
  for(int ieta = 0; ieta < statsNetSum.size(); ieta++) {
    double dEtaAcc = 2. * proton_acceptance.dEta * (ieta + 1);
    file << setw(15) << dEtaAcc << " ";
    double eta_cut = dEtaAcc / 2.;
    file << setw(15) << eta_cut << " ";
    double k2net = statsNetSum[ieta].GetJointCumulantRatio(2, 0, 0, 1);
    double k2neterror = statsNetSum[ieta].GetJointCumulantRatioError(2, 0, 0, 1);
    file << setw(15) << k2net << " ";
    file << setw(15) << k2neterror << " ";
    double k2sum = statsNetSum[ieta].GetJointCumulantRatio(0, 2, 0, 1);
    double k2sumerror = statsNetSum[ieta].GetJointCumulantRatioError(0, 2, 0, 1);
    file << setw(15) << k2sum << " ";
    file << setw(15) << k2sumerror << " ";
    double c2net1 = statsNetSum[ieta].GetJointCumulant(2, 0);
    double c2net2 = statsNetSum[ieta].GetJointCumulant(0, 1);
    double c2err1 = statsNetSum[ieta].GetJointCumulantError(2, 0);
    double c2err2 = statsNetSum[ieta].GetJointCumulantError(0, 1);
    double c2errcov = statsNetSum[ieta].GetJointCumulantsCovariance(2, 0, 0, 1);
    double c2net = (c2net1/c2net2 - 1.) / c2net2;
    double d1 = 1. / c2net2 / c2net2;
    double d2 = -2. * c2net1 / c2net2 / c2net2 / c2net2 + 1. / c2net2 / c2net2;
    double c2neterr = sqrt(d1 * d1 * c2err1 * c2err1 + d2 * d2 * c2err2 * c2err2 + 2. * d1 * d2 * c2errcov);\
    file << setw(15) << c2net << " ";
    file << setw(15) << c2neterr << " ";
    file << endl;
  }
}

// Write D-measure results to file
void WriteProtonResults(const string& filename, vector<SampleMoments::TwoNumberStatistics>& statsNetSum, const ProtonAcceptance& proton_acceptance) {
  ofstream file(filename);
  file << "# Proton fluctuations results" << endl;
  file << "# Events: " << statsNetSum.back().GetNumberOfObservations() << endl;
  WriteProtonResultsSingle(file, statsNetSum, proton_acceptance);
  file.close();
}

// Write D-measure results to file
void WriteProtonResultsAll(const string& filename, vector<vector<SampleMoments::TwoNumberStatistics>>& allstatsNetSum, const vector<ProtonAcceptance>& proton_acceptances) {
  ofstream file(filename);
  file << "# Proton fluctuations results" << endl;
  file << "# Events: " << allstatsNetSum.back().back().GetNumberOfObservations() << endl;
  for (int i = 0; i < proton_acceptances.size(); i++) {
    WriteProtonResultsSingle(file, allstatsNetSum[i], proton_acceptances[i]);
    file << endl << endl;
  }
  file.close();
}

int main(int argc, char* argv[]) {

  cout << "Running FIST sampler version " << FistSampler_VERSION_MAJOR << "." << FistSampler_VERSION_MINOR << endl << endl;

  FistSamplerParameters run_parameters;


  string fileinput = std::string(FistSampler_INPUT_FOLDER) + "/../tasks/Dmeasure/input/input.BCE.ALICE.PbPb.5020.C0-5.EVHRG.ALICEBW";
  fileinput = std::string(FistSampler_INPUT_FOLDER) + "/../tasks/Dmeasure/input/input.BCE.Vc3.ALICE.PbPb.5020.C0-5.ALICEBW";
  fileinput = std::string(FistSampler_INPUT_FOLDER) + "/../tasks/Dmeasure/input/input.BQS.Vc3.ALICE.PbPb.5020.C0-5.ALICEBW";
  //fileinput = std::string(FistSampler_INPUT_FOLDER) + "/../tasks/Dmeasure/input/input.BCE.ALICE.PbPb.5020.C0-5.ALICEBW";

  if (argc > 1) {
    fileinput = string(argv[1]);
    // run_parameters.ReadParametersFromFile(fileinput);
  }

  run_parameters.ReadParametersFromFile(fileinput);
  //run_parameters.hypersurface_file = std::string(FistSampler_INPUT_FOLDER) + "/hydro/AuAu7.7/C70-80/surface_eps_0.26.dat";

  //run_parameters.nevents = 1000;
  {
    int ind = fileinput.size() - 1;
    while (ind >= 0 && fileinput[ind] != '/' && fileinput[ind] != '\\') {
      ind--;
    }
    run_parameters.output_file = fileinput.substr(ind + 1);
    if (run_parameters.output_file.substr(0, 5) == "input") {
      run_parameters.output_file = run_parameters.output_file.substr(6);
    }
    run_parameters.output_file = "ProtonFluctuations." + run_parameters.output_file + ".dat";
  }

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
  if (run_parameters.nevents > 0) {
    cout << "Sampling " << run_parameters.nevents << " events..." << "\n";
  } else {
    cout << "Sampling events indefinitely..." << "\n";
  }

  // Prepare statistics
  // Nch in full space
  SampleMoments::NumberStatistics statsNchtot;

  // Net proton in acceptance regions
  vector<vector<SampleMoments::TwoNumberStatistics>> allstatsNetSum(proton_acceptances.size());
  for(int iacc = 0; iacc < proton_acceptances.size(); iacc++) {
    allstatsNetSum[iacc].resize(proton_acceptances[iacc].eta_iters);
  }

  // Same but for primordial hadrons
  SampleMoments::NumberStatistics statsNchtotPrim;
  vector<vector<SampleMoments::TwoNumberStatistics>> allstatsNetSumPrim(proton_acceptances.size());
  for(int iacc = 0; iacc < proton_acceptances.size(); iacc++) {
    allstatsNetSumPrim[iacc].resize(proton_acceptances[iacc].eta_iters);
  }

  // Loop through the events
  for (long long event_number = 0; event_number < run_parameters.nevents || run_parameters.nevents < 0; ++event_number) {
    // Sample the primordial hadrons
    SimpleEvent evt = evtgen->GetEvent(false);

    // Gather statistics for primordial hadrons
    ProcessEventProtons(evt, allstatsNetSumPrim);

    // Perform the decays, if necessary
    if (lround(run_parameters.parameters["decays"]) != 0) {
      evt = EventGeneratorBase::PerformDecays(evt, TPS);
    }

    // Write the event to file
    if (event_writer != NULL) {
      event_writer->WriteEvent(evt);
    }

    // Process event for the D-measure
    ProcessEventProtons(evt, allstatsNetSum);

    // Periodically print the number of processed events on screen (every 1% or every 1000 events)
    if (run_parameters.nevents < 100 || (event_number + 1) % (run_parameters.nevents / 100) == 0 || (event_number + 1) % 1000 == 0) {
      cout << (event_number + 1) << " ";
      cout.flush();

//      for(int iacc = 0; iacc < nch_acceptances.size(); iacc++) {
//        string filename = "nchacc." + to_string(iacc) + "." + run_parameters.output_file;
//        WriteDmeasureResults(filename, statsNchtot, allstatsNetSum[iacc], nch_acceptances[iacc]);
//      }
      WriteProtonResultsAll("Primordial." + run_parameters.output_file, allstatsNetSumPrim, proton_acceptances);
      WriteProtonResultsAll(run_parameters.output_file, allstatsNetSum, proton_acceptances);
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
