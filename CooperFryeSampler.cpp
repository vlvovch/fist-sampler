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
#include "include/CooperFryeSamplerParameters.h"
#include "include/CooperFryeSamplerHelperFunctions.h"

#include "CooperFryeSamplerConfig.h"

using namespace std;

#ifdef ThermalFIST_USENAMESPACE
using namespace thermalfist;
#endif

using namespace CooperFryeSampler;


int main(int argc, char* argv[]) {

  CooperFryeSamplerParameters run_parameters;

  if (argc > 1) {
    string fileinput = string(argv[1]);
    run_parameters.ReadParametersFromFile(fileinput);
  }

  //run_parameters.hypersurface_file = std::string(CooperFryeSampler_INPUT_FOLDER) + "/hydro/AuAu7.7/C70-80/surface_eps_0.26.dat";


  if (argc > 2) {
    run_parameters.output_file = string(argv[2]);
  }

  // Output the values of all the parameters used
  run_parameters.OutputParameters();

  // Set the random seed
  RandomGenerators::SetSeed(run_parameters.randomseed);

  // Read the Cooper-Frye hypersurface from file
  ParticlizationHypersurface hypersurface;
  ReadHypersurfaceFromFile(run_parameters, hypersurface);

  // Check if the hypersurface is non-empty
  if (hypersurface.size() == 0) {
    std::cout << "Empty hypersurface! Aborting..." << "\n";
    exit(1);
  }

  cout << "Initializing event generator..." << "\n";
  HypersurfaceEventGenerator* evtgen = CreateEventGenerator(run_parameters, hypersurface);
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
