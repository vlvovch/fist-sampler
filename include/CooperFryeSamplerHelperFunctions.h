#ifndef URQMD_PARTICLE_LIST_H
#define URQMD_PARTICLE_LIST_H

#include <set>

#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalParticleSystem.h"
#include "HRGBase/ParticleDecay.h"
#include "CooperFryeSamplerParameters.h"
#include "HypersurfaceReader.h"
#include "HRGEventGenerator/HypersurfaceSampler.h"

namespace CooperFryeSampler {

  // List of PDG codes recognized by UrQMD, all unrecognized will be decayed
  static std::set<long long> UrQMDcodelist =
  {
    // nucleons
    2112, 2212,
    // N*
    12112, 12212, 1214, 2124,
    22112, 22212, 32112, 32212, 2116, 2216,
    12116, 12216, 21214, 22124, 42112, 42212,
    31214, 32124, 1218,  2128,
    // Delta
    1114,  2114,  2214,  2224,
    31114, 32114, 32214, 32224,
    1112,  1212,  2122,  2222,
    11114, 12114, 12214, 12224,
    11112, 11212, 12122, 12222,
    1116,  1216,  2126,  2226,
    21112, 21212, 22122, 22222,
    21114, 22114, 22214, 22224,
    11116, 11216, 12126, 12226,
    1118,  2118,  2218,  2228,
    // Lambda
    3122, 13122, 3124, 23122,
    33122, 13124, 43122, 53122,
    3126, 13126, 23124, 3128, 23126,
    // Sigma
    3112, 3212, 3222,
    3114, 3214, 3224,
    13112, 13212, 13222,
    13114, 13214, 13224,
    23112, 23212, 23222,
    3116,  3216,  3226,
    13116, 13216, 13226,
    23114, 23214, 23224,
    3118,  3218,  3228,
    // Xi
    3312, 3322, 3314, 3324,
    203312, 203322, 13314, 13324,
    103316, 103326, 203316, 203326,
    // Omega
    3334,
    // gamma
    22,
    // mesons
    -211, 111, 211,
    221,223,
    -213, 113, 213,
    9010221,
    311, 321,
    331,
    313, 323,
    333,
    10311, 10321,
    -9000211, 9000111, 9000211,
    10221,
    10313, 10323,
    -20213, 20113, 20213,
    20223,
    40223,
    315, 325,
    -215, 115, 215,
    225, 335,
    20313, 20323,
    -10213, 10113, 10213,
    10223,
    100313, 100323,
    -100213, 100113, 100213,
    100223,
    100333,
    30313, 30323,
    -30213, 30113, 30213,
    30223,
    337,
    421, 411,
    10421, 10411,
    443,
    100443,
    10441,
    431,
    433
  };

  static void ReadHypersurfaceFromFile(CooperFryeSamplerParameters& params, thermalfist::ParticlizationHypersurface& hypersurface)
  {
    int hypersurface_filetype = round(params.parameters["hypersurface_filetype"]);

    if (hypersurface_filetype < 0 || hypersurface_filetype > 2) {
      std::cout << "Invalid hypersurface input file type! Aborting..." << "\n";
      exit(1);
    }

    if (hypersurface_filetype < 2) {
      ReadParticlizationHypersurface(params.hypersurface_file, hypersurface, hypersurface_filetype);
    }
    else if (hypersurface_filetype == 2) {
      ReadParticlizationHypersurfaceMUSIC(params.hypersurface_file, hypersurface);
    }
  }


  static void SetParticleFlags(thermalfist::ThermalParticleSystem* TPS, int decays) {
    for (int ip = 0; ip < TPS->ComponentsNumber(); ++ip) {
      auto& part = TPS->Particle(ip);

      // Boltzmann statistics throughout
      part.UseStatistics(false);

      if (decays == 2) {
        // Strong only
        if (part.DecayType() == thermalfist::ParticleDecayType::Strong)
          part.SetStable(false);
        else
          part.SetStable(true);
      }
      else if (decays == 3) {
        // Strong and electromagnetic
        if (part.DecayType() == thermalfist::ParticleDecayType::Strong
          || part.DecayType() == thermalfist::ParticleDecayType::Electromagnetic)
          part.SetStable(false);
        else
          part.SetStable(true);
      }
      else if (decays == 4) {
        // Strong, electromagnetic, and weak
        if (part.DecayType() == thermalfist::ParticleDecayType::Strong
          || part.DecayType() == thermalfist::ParticleDecayType::Electromagnetic
          || part.DecayType() == thermalfist::ParticleDecayType::Weak)
          part.SetStable(false);
        else
          part.SetStable(true);
      }
      else if (decays == 10) {
        // Set to unstable only if not recognized by UrQMD
        part.SetStable(UrQMDcodelist.count(abs(part.PdgId())));
      }

      // In any case, do not decay charged pions and kaons
      if (part.PdgId() == 211 || part.PdgId() == -211 || part.PdgId() == 321 || part.PdgId() == -321)
        part.SetStable(true);
    }
  }

  thermalfist::HypersurfaceEventGenerator* CreateEventGenerator(CooperFryeSamplerParameters& run_parameters, const thermalfist::ParticlizationHypersurface& hypersurface) {

    std::string list_file = std::string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list.dat";
    if (run_parameters.particle_list_file != "")
      list_file = run_parameters.particle_list_file;
    std::string decays_file = "";
    if (run_parameters.decays_list_file != "")
      decays_file = run_parameters.decays_list_file;

    thermalfist::ThermalParticleSystem* TPS = new thermalfist::ThermalParticleSystem(list_file, decays_file);

    int decays = lround(run_parameters.parameters["decays"]);
    SetParticleFlags(TPS, decays);

    int finite_widths = lround(run_parameters.parameters["finite_widths"]);
    if (finite_widths)
      TPS->SetResonanceWidthIntegrationType(thermalfist::ThermalParticle::eBWconstBR);

    // Create the event generator
    thermalfist::EventGeneratorConfiguration configMC;
    configMC.fEnsemble = thermalfist::EventGeneratorConfiguration::CE;
    configMC.CanonicalB = (lround(run_parameters.parameters["Bcanonical"]) != 0);;
    configMC.CanonicalQ = (lround(run_parameters.parameters["Qcanonical"]) != 0);
    configMC.CanonicalS = (lround(run_parameters.parameters["Scanonical"]) != 0);
    configMC.CanonicalC = (lround(run_parameters.parameters["Ccanonical"]) != 0);

    if (run_parameters.parameters.count("gammaq")) {
      configMC.CFOParameters.gammaq = run_parameters.parameters["gammaq"];
    }

    if (run_parameters.parameters.count("gammaS")) {
      configMC.CFOParameters.gammaS = run_parameters.parameters["gammaS"];
    }

    if (run_parameters.parameters.count("gammaC")) {
      configMC.CFOParameters.gammaC = run_parameters.parameters["gammaC"];
    }

    thermalfist::HypersurfaceEventGenerator* evtgen = new thermalfist::HypersurfaceEventGenerator(
      TPS, configMC, &hypersurface);
    return evtgen;
  }

}

#endif