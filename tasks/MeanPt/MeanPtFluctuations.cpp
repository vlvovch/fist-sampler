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


int main(int argc, char* argv[]) {

  cout << "Running FIST sampler version " << FistSampler_VERSION_MAJOR << "." << FistSampler_VERSION_MINOR << endl << endl;

  FistSamplerParameters run_parameters;


  string fileinput = std::string(FistSampler_INPUT_FOLDER) + "/../tasks/MeanPt/input/input.STAR.AuAu.3GeV.C0-5";

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
  EventGeneratorBase* evtgen2;

  if (fist_sampler_mode == 0 || fist_sampler_mode == 1) {
    evtgen = CreateEventGeneratorFromHypersurface(run_parameters, hypersurface);
  }
  else if (fist_sampler_mode == 2) {
    evtgen = CreateBlastWaveEventGenerator(run_parameters);
    //run_parameters.parameters["BW_betaS"] = 0.80;
    //evtgen2 = CreateBlastWaveEventGenerator(run_parameters);

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

  // Prepare statistics
  SampleMoments::NumberStatistics statsp;
  SampleMoments::NumberStatistics statsMeanPt;
  SampleMoments::NumberStatistics statsCk;

  vector<double> meanPt_entries;
  vector<double> meanPtiPtj_entries;

  // Loop through the events
  for (long long event_number = 0; event_number < run_parameters.nevents; ++event_number) {
    // Sample the primordial hadrons
    SimpleEvent evt;
    //if (event_number < run_parameters.nevents / 2)
    evt = evtgen->GetEvent();
    //else
    //evt = evtgen2->GetEvent();

    // Perform the decays, if necessary
    if (lround(run_parameters.parameters["decays"]) != 0) {
      evt = EventGeneratorBase::PerformDecays(evt, TPS);
    }

    // Write the event to file
    if (event_writer != NULL) {
      //event_writer->WriteEvent(evt);
    }

    // Loop over all the particles
    int Np = 0;
    int Nch = 0;
    double totPt = 0.;
    vector<double> pTi;
    for (auto& part: evt.Particles) {
      // Count protons
      if (part.PDGID == 2212 && part.GetPt() > 0.4 && part.GetPt() < 2.0 && abs(part.GetY()) < 0.5)
        Np++;

      // Count charged particles
      if (evtgen->ThermalModel()->TPS()->ParticleByPDG(part.PDGID).ElectricCharge() != 0) {

        // Momentum cuts
        double pTmin = 0.150, pTmax = 2.000;
        double etamax = 1.0;
        if (part.GetPt() > pTmin && part.GetPt() < pTmax && abs(part.GetEta()) < etamax) {
          Nch++;
          double tPt = part.GetPt();
          totPt += tPt;
          pTi.push_back(tPt);
        }
      }
    }


    double meanPt = totPt / Nch;
    meanPt_entries.push_back(meanPt);

    double ptiptj = 0.;
    for(int i = 0; i < Nch; ++i) {
      for(int j = 0; j < Nch; ++j) {
        if (i != j)
          ptiptj += pTi[i] * pTi[j];
      }
    }
    ptiptj /= Nch * (Nch - 1.);

    meanPtiPtj_entries.push_back(ptiptj);

    statsp.AddObservation(Np);
    statsMeanPt.AddObservation(meanPt);

    // Periodically print the number of processed events on screen
    if (run_parameters.nevents < 100 || (event_number + 1) % (run_parameters.nevents / 100) == 0 || (event_number + 1) % 1000 == 0) {
      cout << (event_number + 1) << " ";
      cout.flush();
    }
  }
  cout << endl;


  // Gather mean pT statistics
  for (long long event_number = 0; event_number < run_parameters.nevents; ++event_number) {
    double Ck = meanPtiPtj_entries[event_number];
    Ck += -2. * statsMeanPt.GetMean() * meanPt_entries[event_number];
    Ck += statsMeanPt.GetMean() * statsMeanPt.GetMean();
    statsCk.AddObservation(Ck);
  }

  // Cleanup
  delete evtgen;
  delete TPS;

  cout << setw(40) << "k1 = " << statsp.GetMean() << " +- " << statsp.GetMeanError() << endl;

  cout << setw(40) << "k2/k1 = " << statsp.GetScaledVariance() << " +- " << statsp.GetScaledVarianceError() << endl;

  cout << setw(40) << "<<p_T>>[GeV] = " << statsMeanPt.GetMean() << " +- " << statsMeanPt.GetMeanError() << endl;

  cout << setw(40) << "<<Delta p_T^2>>[GeV^2] = " << statsMeanPt.GetVariance() << " +- " << statsMeanPt.GetVarianceError() << endl;

  cout << setw(40) << "<<Delta p_T^2>>[GeV] = " << sqrt(statsMeanPt.GetVariance()) << " +- " <<
       statsMeanPt.GetVarianceError()/2./sqrt(statsMeanPt.GetVariance()) << endl;

  cout << setw(40) << "<<Delta p_T_i,j>>[GeV^2] = " << statsCk.GetMean() << " +- " << statsCk.GetMeanError() << endl;

  cout << setw(40) << "sqrt(<<Delta p_T^2>>)/<p_T> = " << sqrt(statsMeanPt.GetVariance()) / statsMeanPt.GetMean() << " +- "
       << statsMeanPt.GetVarianceError() / 2. / sqrt(statsMeanPt.GetVariance()) / statsMeanPt.GetMean() << endl;

  cout << setw(40) << "sqrt(<<Delta p_T_i,j>>)/<p_T> = " << sqrt(abs(statsCk.GetMean())) / statsMeanPt.GetMean() << " +- "
       << statsCk.GetMeanError() / 2. / sqrt(abs(statsCk.GetMean())) / statsMeanPt.GetMean() << endl;


  // Time performace
  double wt2 = get_wall_time();

  cout << "Time per single event: " << (wt2 - wt1) / run_parameters.nevents * 1.e3 << " ms" << "\n";

  return 0;
}
