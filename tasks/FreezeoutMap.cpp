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

#include "../include/HypersurfaceReader.h"

using namespace std;

#ifdef ThermalFIST_USENAMESPACE
using namespace thermalfist;
#endif

using namespace CooperFryeSampler;

double Tmin = 0.0, Tmax = 0.200;
int Titers = 100;
double muBmin = 0.0, muBmax = 0.900;
int muBiters = 90;

string hypersurface_file = "freezeout.dat";
string output_file = "TmuB.7.7.dat";

struct bin {
  double T, muB;
  double weight;
};

int main(int argc, char* argv[]) {

  hypersurface_file =
    "C:\\Users\\volodya\\Documents\\Programming\\ResearchWork\\thirdparty\\ChunHydro\\RHIC7_7\\hydro_results_C0-5\\surface_eps_0.26.dat";
  //hypersurface_file =
  //  "C:\\Users\\volodya\\Documents\\Programming\\ResearchWork\\thirdparty\\ChunHydro\\RHIC7_7\\hydro_results_C70-80\\surface_eps_0.26.dat";

  if (argc > 1) {
    hypersurface_file = string(argv[1]);
  }

  if (argc > 2) {
    output_file = string(argv[2]);
  }

  ParticlizationHypersurface hypersurface;

  ReadParticlizationHypersurfaceMUSIC(hypersurface_file, hypersurface);

  if (!hypersurface.size())
    exit(1);

  vector<vector<bin>> bins;
  double dT = (Tmax - Tmin) / Titers;
  double dmuB = (muBmax - muBmin) / muBiters;
  for (int i = 0; i < Titers; ++i) {
    double T = Tmin + (0.5 + i) * dT;
    bins.push_back(vector<bin>());
    for (int j = 0; j < muBiters; ++j) {
      double muB = Tmin + (0.5 + j) * dmuB;
      bins[i].push_back({ T, muB, 0. });
    }
  }

  double totEnergy = 0.;
  for (const auto& elem : hypersurface) {
    double dVeff = 0.;
    for (int mu = 0; mu < 4; ++mu)
      dVeff += elem.dsigma[mu] * elem.u[mu];

    if (dVeff < 0.)
      continue;

    double En = elem.edens * dVeff;
    totEnergy += En;

    int iT = (int)((elem.T - Tmin) / dT);
    int imuB = (int)((elem.muB - muBmin) / dmuB);
    bins[iT][imuB].weight += En;
  }

  std::cout << "Total energy = " << totEnergy << " GeV" << endl;

  for (auto& el1 : bins) {
    for (auto& el : el1) {
      el.weight /= totEnergy * dT * dmuB;
    }
  }

  ofstream fout(output_file);
  fout << setw(19) << "T[GeV]" << " ";
  fout << setw(19) << "muB[GeV]" << " ";
  fout << setw(19) << "w" << endl;

  for (auto& el1 : bins) {
    for (auto& el : el1) {
      fout << setw(19) << el.T << " ";
      fout << setw(19) << el.muB << " ";
      fout << setw(19) << el.weight << endl;
    }
  }


  return 0;
}
