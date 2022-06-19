#ifndef COOPERFRYE_SAMPLER_PARAMETERS_H
#define COOPERFRYE_SAMPLER_PARAMETERS_H

#include <fstream>
#include <string>
#include <functional>
#include "HRGEventGenerator.h"

#ifdef __linux__
#include <sys/stat.h>
#endif

namespace CooperFryeSampler {

  struct CooperFryeSamplerParameters {
    // Number of events to generate
    long long nevents;

    // Random seed
    unsigned int randomseed;

    // Output file name for the generated events
    std::string output_file;

    // Filename with the Cooper-Frye hypersurface
    std::string hypersurface_file;

    // Filename with the particle list; if blank, the default list is used
    std::string particle_list_file;

    // Filename with the decay channels list; if blank, the decays.dat file in the same folder from the list is used
    std::string decays_list_file;

    // Numeric parameters
    std::map<std::string, double> parameters;

    // Default constructor with the default values
    CooperFryeSamplerParameters(const std::string& input_file = "") :
      nevents(1000),
      randomseed(1),
      output_file("AuAu.7.7.0005.events.dat"),
      hypersurface_file("surface_eps_0.26.dat"),
      particle_list_file(""),
      decays_list_file(""),
      parameters({
        {"Bcanonical",           1},   // the global baryon number is grand-canonical (0) or canonical (1) 
        {"Qcanonical",           0},   // the global electric charge is grand-canonical (0) or canonical (1) 
        {"Scanonical",           0},   // the global strangeness is grand-canonical (0) or canonical (1) 
        {"Ccanonical",           0},   // the global charm is grand-canonical (0) or canonical (1) 
        {"finite_widths",        0},   // Resonance widths treatment: 0 - zero width, 1 - energy-dependent Breit-Wigner
        {"hypersurface_filetype",2},   // 0 - Native ascii, 1 - Native binary, 2 - MUSIC binary
        {"decays", 1},                 // 0 - no decays, 1 - use the stability flags from the list, 2 - strong decays, 3 - strong + electromagnetic, 4 - strong + electromagnetic + weak decays (pi-K-p in the final state), 10 - match the particle list from UrQMD
        {"output_format", 0},          // 0 - native ascii, 1 - tailored for UrQMD afterburner at https://github.com/jbernhard/urqmd-afterburner
        {"b",     0.},                 // Excluded volume parameter for baryons (in fm^3)
        {"radB", -1.},                 // Baryon hard-core radius (in fm). If negative (default), its value is inferred from excluded-volume parameter b
        {"rescaleTmu", 0},             // Rescale the values of T and mu to match energy and baryon densities from hydro, most relevant when EV-HRG model used, less so for Id-HRG. Note that the hypersurface MUST correspond to constant energy density of edens
        {"edens",   0.26},             // The energy density corresponding to the Cooper-Frye hypersurface
        {"use_idealHRG_for_means", 0}, // Use the ideal HRG model when evaluating mean hadron yields, faster initialization at moderate accuracy cost
        {"EVfastmode", 1}              // Use (or not) the fast mode when checking the hard-core overlap of particles. If on, keeps sampling the given particle until no overlap with other particles achieved, instead of rejecting all sampled particles and starting over.
        })
    {
      if (input_file != "")
        ReadParametersFromFile(input_file);
    }

    void ReadParametersFromFile(const std::string& filename) {
      std::ifstream fin(filename);
      if (fin.is_open()) {
        std::cout << "Reading input parameters from file " << filename << "\n";
        std::string var;
        double val;

        while (fin >> var) {
          if (var.size() == 0 || var[0] == '#') {
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
          }

          std::cout << "Reading input parameter " << var << " = ";

          if (var == "output_file") {
            fin >> output_file;
            std::cout << output_file << std::endl;
            continue;
          }

          if (var == "randomseed") {
            fin >> randomseed;
            std::cout << randomseed << std::endl;
            continue;
          }

          if (var == "hypersurface_file") {
            fin >> hypersurface_file;
            std::cout << hypersurface_file << std::endl;
            continue;
          }

          if (var == "particle_list_file") {
            fin >> particle_list_file;
            std::cout << particle_list_file << std::endl;
            continue;
          }

          if (var == "decays_list_file") {
            fin >> decays_list_file;
            std::cout << decays_list_file << std::endl;
            continue;
          }

          fin >> val;
          std::cout << val << std::endl;

          if (var == "nevents") {
            nevents = round(val);
          }
          else {
            parameters[var] = val;
          }
        }
      }
      else {
        std::cout << "Cannot open parameters file!" << "\n";
      }
      std::cout.flush();
    }

    void OutputParameters() {
      const int tabsize = 25;

      std::cout << "Cooper-Frye sampler parameter list:" << "\n";
      std::cout << std::setw(tabsize) << "nevents" << " = " << nevents << "\n";
      std::cout << std::setw(tabsize) << "randomseed" << " = " << randomseed << "\n";

      std::cout << std::setw(tabsize) << "hypersurface_file" << " = " << hypersurface_file << "\n";
      std::cout << std::setw(tabsize) << "particle_list_file" << " = " << particle_list_file << "\n";
      std::cout << std::setw(tabsize) << "decays_list_file" << " = " << decays_list_file << "\n";

      for (auto& el : parameters) {
        std::cout << std::setw(tabsize) << el.first << " = " << el.second << "\n";
      }

      std::cout << std::setw(tabsize) << "output_file" << " = " << output_file << "\n";

      std::cout << std::endl;
    }
  };
}

#endif