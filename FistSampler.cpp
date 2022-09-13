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

#include "include/HypersurfaceReader.h"
#include "include/FistSamplerParameters.h"
#include "include/FistSamplerHelperFunctions.h"

#include "FistSamplerConfig.h"


using namespace std;

#ifdef ThermalFIST_USENAMESPACE
using namespace thermalfist;
#endif

using namespace FistSampler;


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

  ofstream fout_events(run_parameters.output_file);

  // Measure time
  double wt1 = get_wall_time();

  cout << "\n";
  cout << "Sampling " << run_parameters.nevents << " events..." << endl;


  // Loop through the events
  for (long long event_number = 0; event_number < run_parameters.nevents; ++event_number) {
    // Sample the primordial hadrons
    SimpleEvent evt = evtgen->GetEvent();

    // Perform the decays, if necessary
    if (lround(run_parameters.parameters["decays"]) != 0) {
      evt = EventGeneratorBase::PerformDecays(evt, TPS);
    }

    // Write the event to file
    if (event_writer != NULL) {
      event_writer->WriteEvent(evt);
    }

    // Periodically print the number of processed events on screen
    if (run_parameters.nevents < 100 || (event_number + 1) % (run_parameters.nevents / 100) == 0 || (event_number + 1) % 1000 == 0) {
      cout << (event_number + 1) << " ";
      cout.flush();
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
