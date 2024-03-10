#ifndef FIST_SAMPLER_HYPERSURFACE_READER_H
#define FIST_SAMPLER_HYPERSURFACE_READER_H

#include <fstream>
#include <string>
#include <functional>
#include "HRGEventGenerator.h"
#include "FistSamplerConfig.h"

#if defined(__linux__) || defined(__APPLE__)
#include <sys/stat.h>
#endif

namespace FistSampler {

  static void ReadParticlizationHypersurface(const std::string& hypersurface_filename, thermalfist::ParticlizationHypersurface& hypersurface, int binary_file_mode = 0) {
    if (binary_file_mode) {
      // get file size
      std::function<long(std::string)> GetFileSize = [](std::string filename) {
        struct stat stat_buf;
        int rc = stat(filename.c_str(), &stat_buf);
        return rc == 0 ? stat_buf.st_size : -1;
      };


      const size_t filesize = GetFileSize(hypersurface_filename);

      std::ifstream fin(hypersurface_filename, std::ios::in | std::ios::binary);

      if (fin.is_open()) {
        const size_t count = filesize / sizeof(thermalfist::ParticlizationHypersurfaceElement);


        hypersurface = thermalfist::ParticlizationHypersurface(count);
        fin.read(reinterpret_cast<char*>(&hypersurface[0]), count * sizeof(thermalfist::ParticlizationHypersurfaceElement));
      }

      fin.close();
    }
    else {
      std::ifstream fin(hypersurface_filename);
      if (fin.is_open()) {

        hypersurface.clear();

        double tau, x, y, eta;
        double dsigma[4];
        double u[4];
        double T, muB, muQ, muS;
        int Nelem = 0;
        while (fin >> tau >> x >> y >> eta) {
          if (Nelem % 100000 == 0)
            std::cout << Nelem << " ";
          for (int ii = 0; ii < 4; ++ii)
            fin >> dsigma[ii];
          for (int ii = 0; ii < 4; ++ii)
            fin >> u[ii];
          fin >> T >> muB >> muQ >> muS;
          double tmp;
          for (int ii = 0; ii < 11; ++ii)
            fin >> tmp;

          Nelem++;

          hypersurface.push_back({
            tau, x, y, eta, 
            dsigma[0], dsigma[1], dsigma[2], dsigma[3],
            u[0], u[1], u[2], u[3],
            T, muB, muQ, muS
            });
        }
      }

      fin.close();
    }
  }

  static void WriteHypersurfaceToBinaryFile(const thermalfist::ParticlizationHypersurface& hypersurface, const std::string& filename) {
    std::ofstream fout(filename, std::ios::out | std::ios::binary);
    fout.write(reinterpret_cast<const char*>(&hypersurface[0]), hypersurface.size() * sizeof(thermalfist::ParticlizationHypersurfaceElement));
    fout.close();
  }

  static void ConvertFreezeoutHypersurface(const std::string& hypersurface_filename, const std::string& hypersurface_filename_binary) {
    thermalfist::ParticlizationHypersurface hypersurface;
    ReadParticlizationHypersurface(hypersurface_filename, hypersurface, 0);
    std::cout << hypersurface.size() << " ";
    WriteHypersurfaceToBinaryFile(hypersurface, hypersurface_filename_binary);
  }

  // Following the iSS sampler https://github.com/chunshen1987/iSS/blob/main/src/readindata.cpp#L626
  static void ReadParticlizationHypersurfaceMUSIC(const std::string& hypersurface_filename, 
    thermalfist::ParticlizationHypersurface& hypersurface
    ) {
    std::ifstream fin(hypersurface_filename, std::ios::binary);

    if (!fin.is_open()) {
      std::cout << "Hypersurface file " << hypersurface_filename << " not found!" << std::endl;
      std::string hypersurface_filename2 = std::string(FistSampler_INPUT_FOLDER) + "/" + hypersurface_filename;
      std::cout << "Trying " << hypersurface_filename2 << std::endl;
      fin.open(hypersurface_filename2, std::ios::binary);
      if (!fin.is_open()) {
        std::cout << "Hypersurface file " << hypersurface_filename2 << " not found!" << std::endl;
        std::cout << "Aborting..." << std::endl;
        return;
      }
    }

    hypersurface.clear();

    while (!fin.eof()) {
      float arr[34];
      for (int i = 0; i < 34; i++) {
        float temp = 0.;
        fin.read(reinterpret_cast<char*>(&temp), sizeof(float));
        arr[i] = temp;
      }

      thermalfist::ParticlizationHypersurfaceElement elem;
      elem.tau  = arr[0];
      elem.x = arr[1];
      elem.y = arr[2];
      elem.eta = arr[3];

      double da0 = arr[4];
      double da1 = arr[5];
      double da2 = arr[6];
      double da3 = arr[7];
      double uu0 = arr[8];
      double uu1 = arr[9];
      double uu2 = arr[10];
      double uu3 = arr[11];

      elem.edens = arr[12] / thermalfist::xMath::GeVtoifm();
      elem.rhoB  = arr[29];

      elem.T   = arr[13] / thermalfist::xMath::GeVtoifm();
      elem.muB = arr[14] / thermalfist::xMath::GeVtoifm();
      elem.muS = arr[15] / thermalfist::xMath::GeVtoifm();
      elem.muQ = arr[16] / thermalfist::xMath::GeVtoifm();

      double dVeff1 = elem.tau * (da0 * uu0 + da1 * uu1 + da2 * uu2 + da3 * uu3 / elem.tau);
      double umod1 = uu0 * uu0 - uu1 * uu1 - uu2 * uu2 - uu3 * uu3;

      double cosheta = std::cosh(elem.eta), sinheta = std::sinh(elem.eta);

      // Transform u and dSigma to Cartesian coordinates
      elem.u[0] = uu0 * cosheta + uu3 * sinheta;
      elem.u[1] = uu1;
      elem.u[2] = uu2;
      elem.u[3] = uu0 * sinheta + uu3 * cosheta;

      elem.dsigma[0] = elem.tau * cosheta * da0 - sinheta * da3;
      elem.dsigma[1] = elem.tau * da1;
      elem.dsigma[2] = elem.tau * da2;
      elem.dsigma[3] = -(elem.tau * sinheta * da0 - cosheta * da3);

      double dVeff2 = 0.;
      for (int mu = 0; mu < 4; ++mu)
        dVeff2 += elem.dsigma[mu] * elem.u[mu];

      double umod2 = elem.u[0] * elem.u[0] - elem.u[1] * elem.u[1] - elem.u[2] * elem.u[2] - elem.u[3] * elem.u[3];

      // Shear stress tensor and pressure
      //float* pi = &arr[18];
      // First, in the Milne coordiantes
      double piMilne00 = arr[18] / thermalfist::xMath::GeVtoifm();
      double piMilne01 = arr[19] / thermalfist::xMath::GeVtoifm();
      double piMilne02 = arr[20] / thermalfist::xMath::GeVtoifm();
      double piMilne03 = arr[21] / thermalfist::xMath::GeVtoifm();
      double piMilne11 = arr[22] / thermalfist::xMath::GeVtoifm();
      double piMilne12 = arr[23] / thermalfist::xMath::GeVtoifm();
      double piMilne13 = arr[24] / thermalfist::xMath::GeVtoifm();
      double piMilne22 = arr[25] / thermalfist::xMath::GeVtoifm();
      double piMilne23 = arr[26] / thermalfist::xMath::GeVtoifm();
      double piMilne33 = arr[27] / thermalfist::xMath::GeVtoifm();
      // Now transfer to Cartesian coordinates and index44
      // Transform to Cartesian coordinates as per https://github.com/chunshen1987/iSS/blob/main/src/iSS.cpp#L248
      // pi[0] == pi_tz[0][0]
      elem.pi[0] = piMilne00 * cosheta * cosheta +
              2. * piMilne03 * cosheta * sinheta +
              piMilne33 * sinheta * sinheta;
      // pi[1] == pi_tz[0][1]
      elem.pi[1] = piMilne01 * cosheta + piMilne13 * sinheta;
      // pi[2] == pi_tz[1][1]
      elem.pi[2] = piMilne11;
      // pi[3] == pi_tz[0][2]
      elem.pi[3] = piMilne02 * cosheta + piMilne23 * sinheta;
      // pi[4] == pi_tz[1][2]
      elem.pi[4] = piMilne12;
      // pi[5] == pi_tz[2][2]
      elem.pi[5] = piMilne22;
      // pi[6] == pi_tz[0][3]
      elem.pi[6] = piMilne00 * cosheta * sinheta +
                   piMilne03 * (cosheta*cosheta + sinheta*sinheta) +
                   piMilne33 * sinheta * cosheta;
      // pi[7] == pi_tz[1][3]
      elem.pi[7] = piMilne01 * sinheta + piMilne13 * cosheta;
      // pi[8] == pi_tz[2][3]
      elem.pi[8] = piMilne02 * sinheta + piMilne23 * cosheta;
      // pi[9] == pi_tz[3][3]
      elem.pi[9] = piMilne00 * sinheta * sinheta +
                   2. * piMilne03 * sinheta * cosheta +
                   piMilne33 * cosheta * cosheta;

      // Now calculate the pressure
      // array[17] == (e + P) / T
      elem.P = (arr[17] * elem.T - elem.edens);

      if (!fin.eof() && elem.T > 0.01) // As in iSS
        hypersurface.push_back(elem);
    }
  }
}

#endif