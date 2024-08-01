#ifndef FIST_SAMPLER_HELPER_FUNCTIONS_H
#define FIST_SAMPLER_HELPER_FUNCTIONS_H

#include <set>

#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalParticleSystem.h"
#include "HRGBase/ParticleDecay.h"
#include "FistSamplerParameters.h"
#include "HypersurfaceReader.h"
#include "HRGEventGenerator/HypersurfaceSampler.h"
#include "HRGEV/ExcludedVolumeHelper.h"
#include "HRGBase/xMath.h"

namespace FistSampler {

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

  static void ReadHypersurfaceFromFile(FistSamplerParameters& params, thermalfist::ParticlizationHypersurface& hypersurface)
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


  thermalfist::HypersurfaceEventGenerator* CreateEventGeneratorFromHypersurface(FistSamplerParameters& run_parameters, const thermalfist::ParticlizationHypersurface& hypersurface) {

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
    configMC.CanonicalB = (lround(run_parameters.parameters["Bcanonical"]) != 0);
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

    // Excluded volume for baryons
    double b = run_parameters.parameters["b"];
    if (b > 0.0 && lround(run_parameters.parameters["use_idealHRG_for_means"]) == 0) {
      configMC.fModelType = thermalfist::EventGeneratorConfiguration::CrosstermsEV;
    }

    configMC.bij = std::vector<std::vector<double> >(TPS->ComponentsNumber(), std::vector<double>(TPS->ComponentsNumber(), 0.));
    for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
      for (int j = 0; j < TPS->ComponentsNumber(); ++j) {
        if (TPS->Particle(i).BaryonCharge() * TPS->Particle(j).BaryonCharge() == 1) {
          configMC.bij[i][j] = b;
        }
        else {
          configMC.bij[i][j] = 0.;
        }
      }
    }

    thermalfist::HypersurfaceEventGenerator* evtgen;
    if (b == 0.0)
      evtgen = new thermalfist::HypersurfaceEventGenerator(
      TPS, configMC, &hypersurface);
    else {
      evtgen = new thermalfist::HypersurfaceEventGeneratorEVHRG(
        TPS, configMC, &hypersurface);
      static_cast<thermalfist::HypersurfaceEventGeneratorEVHRG*>(evtgen)->SetExcludedVolume(b);
      double radB = run_parameters.parameters["radB"];
      if (radB < 0.0)
        radB = thermalfist::CuteHRGHelper::rv(b);
      static_cast<thermalfist::HypersurfaceEventGeneratorEVHRG*>(evtgen)->SetBaryonRadius(radB);
    }

    evtgen->SetRescaleTmu((lround(run_parameters.parameters["rescaleTmu"]) != 0), run_parameters.parameters["edens"]);
    evtgen->SetEVUseSPR((lround(run_parameters.parameters["EVfastmode"]) != 0));

    bool shear_correction = (lround(run_parameters.parameters["shear_correction"]) != 0);
    evtgen->SetShearCorrection(shear_correction);

    return evtgen;
  }

  void CreateSiemensRasmussenHubbleHypersurface(FistSamplerParameters& run_parameters, thermalfist::ParticlizationHypersurface& hypersurface) {
    double T = 0.070; // GeV
    if (run_parameters.parameters.count("T"))
      T = run_parameters.parameters["T"];

    double muB = 0.8721; // GeV
    if (run_parameters.parameters.count("muB"))
      muB = run_parameters.parameters["muB"];

    double muQ = -0.0211; // GeV
    if (run_parameters.parameters.count("muQ"))
      muQ = run_parameters.parameters["muQ"];

    double muS = 0.1984; // GeV
    if (run_parameters.parameters.count("muS"))
      muS = run_parameters.parameters["muS"];

    double gammaS = 0.052;
    if (run_parameters.parameters.count("gammaS"))
      gammaS = run_parameters.parameters["gammaS"];

    double H = 0.097; // fm^-1
    if (run_parameters.parameters.count("SRH_H"))
      H = run_parameters.parameters["SRH_H"];

    double R = 6.1; // fm
    if (run_parameters.parameters.count("SRH_R"))
      R = run_parameters.parameters["SRH_R"];

    double t = 15.0; // fm (freeze-out is isochronous, so the value does not matter)

    int iterR = 100;
    if (run_parameters.parameters.count("SRH_iterR"))
      iterR = lround(run_parameters.parameters["SRH_iterR"]);

    int iterCosTh = 100;
    if (run_parameters.parameters.count("SRH_iterCosTh"))
      iterCosTh = lround(run_parameters.parameters["SRH_iterCosTh"]);

    int iterPhi = 100;
    if (run_parameters.parameters.count("SRH_iterPhi"))
      iterPhi = lround(run_parameters.parameters["SRH_iterPhi"]);

    double dR = R / iterR;
    double dCosTh = 2. / iterCosTh;
    double dPhi = 2. * thermalfist::xMath::Pi() / iterPhi;

    hypersurface.clear();

    for (int ir = 0; ir < iterR; ++ir) {
      double r = dR * (0.5 + ir);
      for (int icth = 0; icth < iterCosTh; ++icth) {
        double costh = -1. + dCosTh * (0.5 + icth);
        double sinth = sqrt(1. - costh * costh);
        for (int iphi = 0; iphi < iterPhi; ++iphi) {
          double phi = dPhi * (0.5 + iphi);
          double cosphi = cos(phi);
          double sinphi = sin(phi);

          thermalfist::ParticlizationHypersurfaceElement elem;
          elem.T = T;
          elem.muB = muB;
          elem.muQ = muQ;
          elem.muS = muS;

          double x = r * cosphi * sinth;
          double y = r * sinphi * sinth;
          double z = r * costh;
          elem.tau = sqrt(t * t - z * z);
          elem.x = x;
          elem.y = y;
          elem.eta = 0.5 * log((t - z) / (t + z));

          elem.dsigma[0] = r * r * dCosTh * dPhi * dR;
          elem.dsigma[1] = elem.dsigma[2] = elem.dsigma[3] = 0;

          double vr = tanh(H * r);
          double gammav = 1. / sqrt(1. - vr * vr);

          elem.u[0] = gammav;
          elem.u[1] = gammav * vr * cosphi * sinth;
          elem.u[2] = gammav * vr * sinphi * sinth;
          elem.u[3] = gammav * vr * costh;

          elem.rhoB = elem.edens = 0.;

          hypersurface.push_back(elem);
        }
      }
    }
  }

  thermalfist::CylindricalBlastWaveEventGenerator* CreateBlastWaveEventGenerator(FistSamplerParameters& run_parameters) {

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
    configMC.CanonicalB = (lround(run_parameters.parameters["Bcanonical"]) != 0);
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

    // Excluded volume for baryons
    double b = run_parameters.parameters["b"];
    if (b > 0.0 && lround(run_parameters.parameters["use_idealHRG_for_means"]) == 0) {
      configMC.fModelType = thermalfist::EventGeneratorConfiguration::CrosstermsEV;
    }

    configMC.bij = std::vector<std::vector<double> >(TPS->ComponentsNumber(), std::vector<double>(TPS->ComponentsNumber(), 0.));
    for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
      for (int j = 0; j < TPS->ComponentsNumber(); ++j) {
        if (TPS->Particle(i).BaryonCharge() * TPS->Particle(j).BaryonCharge() == 1) {
          configMC.bij[i][j] = b;
        }
        else {
          configMC.bij[i][j] = 0.;
        }
      }
    }

    // Chemical freeze-out parameters
    configMC.CFOParameters.T = 0.160; // GeV
    if (run_parameters.parameters.count("T"))
      configMC.CFOParameters.T = run_parameters.parameters["T"];

    configMC.CFOParameters.muB = 0.000; // GeV
    if (run_parameters.parameters.count("muB"))
      configMC.CFOParameters.muB = run_parameters.parameters["muB"];

    configMC.CFOParameters.muQ = 0.000; // GeV
    if (run_parameters.parameters.count("muQ"))
      configMC.CFOParameters.muQ = run_parameters.parameters["muQ"];

    configMC.CFOParameters.muS = 0.000; // GeV
    if (run_parameters.parameters.count("muS"))
      configMC.CFOParameters.muS = run_parameters.parameters["muS"];

    configMC.CFOParameters.muC = 0.000; // GeV
    if (run_parameters.parameters.count("muC"))
      configMC.CFOParameters.muC = run_parameters.parameters["muC"];

    double dVdy = 4000.; // fm^3
    if (run_parameters.parameters.count("dVdy"))
      dVdy = run_parameters.parameters["dVdy"];

    // Blast-wave parameters
    double BW_Tkin = configMC.CFOParameters.T;
    if (run_parameters.parameters.count("BW_Tkin"))
      BW_Tkin = run_parameters.parameters["BW_Tkin"];

    double BW_betaS = 0.77;
    if (run_parameters.parameters.count("BW_betaS"))
      BW_betaS = run_parameters.parameters["BW_betaS"];

    double BW_n = 0.36;
    if (run_parameters.parameters.count("BW_n"))
      BW_n = run_parameters.parameters["BW_n"];

    double BW_etamax = 4.8;
    if (run_parameters.parameters.count("BW_etamax"))
      BW_etamax = run_parameters.parameters["BW_etamax"];

    double BW_rmax = 9.0;
    if (run_parameters.parameters.count("BW_rmax"))
      BW_rmax = run_parameters.parameters["BW_rmax"];

    double Veff = dVdy * 2. * BW_etamax;
    configMC.CFOParameters.V = configMC.CFOParameters.SVc = Veff;

    // In case we do canonical ensemble, we will calculate total conserved charges grand-canonically assuming the system is large enough
    // or has zero conserved charges (LHC)
    configMC.fUseGCEConservedCharges = true;

    thermalfist::CylindricalBlastWaveEventGenerator* evtgen = new thermalfist::CylindricalBlastWaveEventGenerator(
      TPS, configMC,
      BW_Tkin, BW_betaS, BW_etamax, BW_n, BW_rmax);
    evtgen->SetEVUseSPR((lround(run_parameters.parameters["EVfastmode"]) != 0));
    evtgen->RecalculateTotalConservedNumbers();

    return evtgen;
  }

}

#endif