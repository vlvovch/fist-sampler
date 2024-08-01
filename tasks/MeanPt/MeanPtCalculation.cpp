// Mean pT fluctuations
// Test only on 3 GeV Au-Au collisions

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

// List of species
vector<long long> species_list = {211, 321, 2212, 3122, 3312, 3334};

// Rapidity cut
double ycut = 0.5;


// Process a single event
void ProcessEventForMeanPts(const SimpleEvent& event,
                            SampleMoments::NumberStatistics& statsNchPt,
                            vector<SampleMoments::NumberStatistics>& allstatsPts) {
  int Nchtot = 0;
  double totPtNch = 0.;
  vector<int> Nps(allstatsPts.size(),0);
  vector<double> totPts(allstatsPts.size(),0.);

  // Iterate over all particles in the event
  for (const auto& particle : event.Particles) {
    // Check if the particle is a charged hadron
    int Q = parts_pdg.ParticleByPDG(particle.PDGID).ElectricCharge();
    if (Q != 0) {
      // Check if particle is inside rapidity cut
      if (abs(particle.GetY()) < ycut) {
        totPtNch += particle.GetPt();
        statsNchPt.AddObservation(particle.GetPt());
        Nchtot++;
      }
    }

    for(int ipdg = 0; ipdg < species_list.size(); ipdg++) {
      auto pdg = species_list[ipdg];
      if (particle.PDGID == pdg) {
        totPts[ipdg] += particle.GetPt();
        allstatsPts[ipdg].AddObservation(particle.GetPt());
        Nps[ipdg]++;
      }
    }
  }

  // Add observables to the statistics
  //statsNchPt.AddObservation(totPtNch / Nchtot);
//  for(int ipdg = 0; ipdg < species_list.size(); ipdg++) {
//    allstatsPts[ipdg].AddObservation(totPts[ipdg] / Nps[ipdg]);
//  }
}

void WriteMeanPtResults(ostream& ostr,
        SampleMoments::NumberStatistics& statsNchPt,
        vector<SampleMoments::NumberStatistics>& allstatsPts) {
  ostr << setw(15) << "Entries" << " "
       << setw(15) << "Species" << " "
       << setw(15) << "<pT>[MeV]" << " "
       << setw(15) << "Error" << " ";
  ostr << endl;

  for(int ipdg = 0; ipdg < species_list.size(); ipdg++) {
    ostr << setw(15) << allstatsPts[ipdg].GetNumberOfObservations() << " ";
    ostr << setw(15) << species_list[ipdg] << " ";
    ostr << setw(15) << 1.e3 * allstatsPts[ipdg].GetMean() << " ";
    ostr << setw(15) << 1.e3 * allstatsPts[ipdg].GetMeanError() << " ";
    ostr << endl;
  }

  ostr << setw(15) << statsNchPt.GetNumberOfObservations() << " ";
  ostr << setw(15) << "Charged" << " ";
  ostr << setw(15) << 1.e3 * statsNchPt.GetMean() << " ";
  ostr << setw(15) << 1.e3 * statsNchPt.GetMeanError() << " ";
  ostr << endl;
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
  SampleMoments::NumberStatistics statsNchPt;
  vector<SampleMoments::NumberStatistics> allstatsPts(species_list.size());

  // Loop through the events
  for (long long event_number = 0; event_number < run_parameters.nevents; ++event_number) {
    // Sample the primordial hadrons
    SimpleEvent evt = evtgen->GetEvent(false);

    // Gather statistics for primordial hadrons
    //ProcessEventForMeanPts(evt, statsNchPt, allstatsPts);

    // Perform the decays, if necessary
    if (lround(run_parameters.parameters["decays"]) != 0) {
      evt = EventGeneratorBase::PerformDecays(evt, TPS);
    }

    // Write the event to file
    if (event_writer != NULL) {
      event_writer->WriteEvent(evt);
    }

    // Process event for the mean Pt
    ProcessEventForMeanPts(evt, statsNchPt, allstatsPts);

    // Periodically print the number of processed events on screen (every 1% or every 1000 events)
    if (run_parameters.nevents < 100 || (event_number + 1) % (run_parameters.nevents / 1) == 0 || (event_number + 1) % 1000 == 0) {
      cout << (event_number + 1) << " " << endl;
      cout.flush();

      // Print the statistics
      WriteMeanPtResults(cout, statsNchPt, allstatsPts);
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
