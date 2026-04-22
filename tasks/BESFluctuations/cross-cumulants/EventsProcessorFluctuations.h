#ifndef EVENTS_PROCESSOR_FLUCTUATIONS_H
#define EVENTS_PROCESSOR_FLUCTUATIONS_H

// Gather various statistics from events for fluctuation analysis

#include "sample-moments/include/NumberStatistics.h"
#include "sample-moments/include/TwoNumberStatistics.h"
#include "HRGEventGenerator.h"

using namespace std;
using namespace thermalfist;
using namespace SampleMoments;

struct stval {
  double sum;
  double sumsqr;
  int nev;
  double wsum, w2sum;
  stval() { sum = sumsqr = 0.; nev = 0; wsum = w2sum = 0.; }
  //double getAv() const { return sum / nev; }
  double getAv() const { return sum / wsum; }
  double getSigma() const { return sqrt((sumsqr / wsum - getAv() * getAv()) / (wsum * wsum / w2sum - 1.)); }
  //double getSigma() const { return sqrt((sumsqr - sum*sum/nev)/(nev-1.)/nev); }
  double getSigma2() const { return sqrt(sum) / nev; }
};

struct density
{
  double xst;
  double step;
  int counter;
  //double xend;
  int sz;
  std::vector<stval> data;
  std::vector<double> tsum;
  int valcount;
  int events;
  density()
  {
  }
  density(double xleft, double xright, int numb)
  {
    if (xright < xleft) printf("Right limit is less than left!\n");
    sz = numb;
    counter = 0;
    step = (xright - xleft) / sz;
    xst = xleft + 0.5 * step;
    data.resize(sz);
    tsum.resize(sz);
    for (int i = 0; i < tsum.size(); ++i) tsum[i] = 0.;
    valcount = 0;
    events = 0;
  }
  void insert(double value)
  {
    int tind = (int)((value - (xst - 0.5 * step)) / step);
    if (tind >= 0 && tind < sz)
    {
      tsum[tind]++;
    }
    counter++;
  }
  void updateEvent(double weight = 1.)
  {
    for (int i = 0; i < tsum.size(); ++i) {
      data[i].sum += weight * tsum[i];
      data[i].sumsqr += weight * tsum[i] * tsum[i];
      data[i].nev++;
      data[i].wsum += weight;
      data[i].w2sum += weight * weight;
      tsum[i] = 0.;
    }
    events++;
  }
  double GetEntry(int ind) const {
    if (ind >= 0 && ind < data.size()) {
      return data[ind].getAv() / step;
    }
    return 0.;
  }
  double GetEntryError(int ind) const {
    if (ind >= 0 && ind < data.size()) {
      return data[ind].getSigma() / step;
    }
    return 0.;
  }
  double GetX(int ind) const {
    if (ind >= 0 && ind < data.size()) {
      return xst + ind * step;
    }
    return 0.;
  }
  std::vector<double> GetXVector() const {
    std::vector<double> ret(data.size());
    for (int i = 0; i < data.size(); ++i)
      ret[i] = xst + i * step;
    return ret;
  }
  std::vector<double> GetYVector() const {
    std::vector<double> ret(data.size(), 0.);
    if (events > 0)
      for (int i = 0; i < data.size(); ++i)
        ret[i] = data[i].getAv() / step;
    return ret;
  }
  std::vector<double> GetYErrorVector() const {
    std::vector<double> ret(data.size(), 0.);
    if (events > 0)
      for (int i = 0; i < data.size(); ++i)
        ret[i] = data[i].getSigma() / step;
    return ret;
  }
};


class EventsProcessorFluctuations
{
public:
  int m_nsubs;
  double m_dY;
  long long nevents;

  // Baryon vs eta_s/rapidity
  vector<vector<NumberStatistics>> statsB;
  vector<vector<NumberStatistics>> statsBY; 

  // Baryon/antibaryon correlators
  vector<TwoNumberStatistics> twostatsB, twostatsp, twostatsnetBp;

  // Net-baryon number
  vector<SampleMoments::NumberStatistics> statsnetB;
  // Net/sum baryon number correlators
  vector<SampleMoments::TwoNumberStatistics> twostatsnetsumB;

  // Net-proton vs eta_s/rapidity
  vector<vector<NumberStatistics>> statsp;
  vector<vector<NumberStatistics>> statspY;

  // Light nuclei rapidity acceptance
  vector<vector<NumberStatistics>> statsdY;
  vector<vector<NumberStatistics>> statstY;
  vector<vector<NumberStatistics>> statsHe3Y;
  vector<vector<NumberStatistics>> statsHe4Y;

  // ALICE identified
  vector<vector<NumberStatistics>> statspALICE;
  vector<vector<NumberStatistics>> statsKALICE;
  vector<vector<NumberStatistics>> statspiALICE;
  vector<vector<NumberStatistics>> statsLALICE;
  vector<vector<NumberStatistics>> statschALICE;

  // STAR identified
  vector<vector<NumberStatistics>> statspSTAR;
  vector<SampleMoments::TwoNumberStatistics> twostatspSTAR;
  vector<SampleMoments::TwoNumberStatistics> twostatsnetpnetkSTAR;
  vector<SampleMoments::TwoNumberStatistics> twostatsnetQnetkSTAR;
  vector<SampleMoments::TwoNumberStatistics> twostatsnetQnetpSTAR;
  vector<SampleMoments::TwoNumberStatistics> twostatsnetpinetpSTAR;
  vector<SampleMoments::TwoNumberStatistics> twostatsnetpinetkSTAR;

  // Pointers to all the stats
  vector<vector<vector<NumberStatistics>>*> statsRapidity;
  vector<vector<vector<NumberStatistics>>*> statsALICE;
  vector<vector<vector<NumberStatistics>>*> statsSTAR;


  // A la https://arxiv.org/pdf/2205.10030.pdf
  vector<SampleMoments::TwoNumberStatistics> twostatsBQKitazawa;
  vector<SampleMoments::TwoNumberStatistics> twostatsBQKitazawapT;
  vector<SampleMoments::TwoNumberStatistics> twostatsBQKitazawaEta;
  vector<SampleMoments::TwoNumberStatistics> twostatsBQKitazawaEtaMix;
  vector<SampleMoments::TwoNumberStatistics> twostatsBQKitazawapTEta;
  vector<SampleMoments::TwoNumberStatistics> twostatsBQKitazawapTEtaMix;
  vector<SampleMoments::TwoNumberStatistics> twostatsBQKitazawaEtaS;

  vector<SampleMoments::TwoNumberStatistics> twostatsBQKitazawaPrim;
  vector<SampleMoments::TwoNumberStatistics> twostatsBQKitazawaPrimEta;
  vector<SampleMoments::TwoNumberStatistics> twostatsBQKitazawaPrimEtaMix;

  // STAR net-kaon/net-Lambda
  vector<SampleMoments::TwoNumberStatistics> twostatsnetsumKaon;
  vector<SampleMoments::TwoNumberStatistics> twostatsnetsumLambda;
  vector<SampleMoments::TwoNumberStatistics> twostatsnetsumKaonNoCut;
  vector<SampleMoments::TwoNumberStatistics> twostatsnetsumLambdaNoCut;

  //vector<TwoNumberStatistics> statsNpNm;
  TwoNumberStatistics statsNpNmtot;

  // pT spectra at mid-rapidity
  int m_iterspT;
  double m_dpT;
  vector<NumberStatistics> statsppT, statsKpT, statspipT;
  NumberStatistics statstotppT, statstotKpT, statstotpipT;
  NumberStatistics statspMeanpT, statsKMeanpT, statspiMeanpT;

  vector<vector<NumberStatistics>*> statspTs;

  NumberStatistics statsdNchdEta;
  NumberStatistics statsNpfull;
  NumberStatistics statsNdfull;
  NumberStatistics statsNtfull;

  // Angular correlations of protons
  density protAngular;

  EventsProcessorFluctuations(int nsubs, double dY, int iterspT, double dpT) {
    nevents = 0;

    m_nsubs = nsubs;
    m_dY = dY;
    m_iterspT = iterspT;
    m_dpT = dpT;

    statsRapidity = { 
      &statsB,
      &statsp,
      &statsBY,
      &statspY,
      & statsdY,
      & statstY,
      & statsHe3Y,
      & statsHe4Y
    };

    statsALICE = {
      &statspALICE,
      &statsKALICE, 
      &statspiALICE,
      &statsLALICE, 
      &statschALICE
    };

    statsSTAR = {
      &statspSTAR
    };

    for (auto& stats : statsRapidity) {
      (*stats) = vector<vector<NumberStatistics>>(3, vector<NumberStatistics>(m_nsubs / 2));
    }

    twostatsB     = vector<TwoNumberStatistics>(m_nsubs / 2);
    twostatsp     = vector<TwoNumberStatistics>(m_nsubs / 2);
    twostatsnetBp = vector<TwoNumberStatistics>(m_nsubs / 2);

    statsnetB = vector<SampleMoments::NumberStatistics>(m_nsubs / 2);
    twostatsnetsumB = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);

    twostatspSTAR = vector<SampleMoments::TwoNumberStatistics>(m_nsubs * (10 / 2));
    twostatsnetpnetkSTAR = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsnetQnetkSTAR = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsnetQnetpSTAR = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsnetpinetpSTAR = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsnetpinetkSTAR = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);



    twostatsBQKitazawa      = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsBQKitazawaEta   = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsBQKitazawaEtaMix = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsBQKitazawapT    = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsBQKitazawapTEta = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsBQKitazawapTEtaMix = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsBQKitazawaEtaS  = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsBQKitazawaPrim  = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsBQKitazawaPrimEta = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsBQKitazawaPrimEtaMix = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);


    twostatsnetsumKaon   = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsnetsumLambda = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);

    twostatsnetsumKaonNoCut   = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);
    twostatsnetsumLambdaNoCut = vector<SampleMoments::TwoNumberStatistics>(m_nsubs / 2);

    for (auto& stats : statsALICE) {
      (*stats) = vector<vector<NumberStatistics>>(3, vector<NumberStatistics>(m_nsubs * (10 / 2)));
    }

    for (auto& stats : statsSTAR) {
      (*stats) = vector<vector<NumberStatistics>>(3, vector<NumberStatistics>(m_nsubs * (10 / 2)));
    }

    //statsNpNm = vector<TwoNumberStatistics>(m_nsubs * (10 / 2));
    TwoNumberStatistics statsNpNmtot;

    statspTs = { &statsppT, &statsKpT, &statspipT };
    for (auto& stats : statspTs) {
      (*stats) = vector<NumberStatistics>(m_iterspT);
    }

    statsdNchdEta = NumberStatistics();

    protAngular = density(-xMath::Pi(), xMath::Pi(), 1000);
  }

  void ProcessEvent(const SimpleEvent& evt, ThermalParticleSystem& TPS) {
    nevents++;

    // eta_s
    vector<int> NBps(m_nsubs / 2, 0), NBms(m_nsubs / 2, 0);
    vector<int> Npps(m_nsubs / 2, 0), Npms(m_nsubs / 2, 0);
    // rapidity
    vector<int> NBpsY(m_nsubs / 2, 0), NBmsY(m_nsubs / 2, 0);
    vector<int> NppsY(m_nsubs / 2, 0), NpmsY(m_nsubs / 2, 0);
    vector<int> NdpsY(m_nsubs / 2, 0), NdmsY(m_nsubs / 2, 0);
    vector<int> NtpsY(m_nsubs / 2, 0), NtmsY(m_nsubs / 2, 0);
    vector<int> NHe3psY(m_nsubs / 2, 0), NHe3msY(m_nsubs / 2, 0);
    vector<int> NHe4psY(m_nsubs / 2, 0), NHe4msY(m_nsubs / 2, 0);

    // ALICE identified
    double pmin = 0.6, pmax = 1.5;
    vector<int> NppALICE(m_nsubs * (10 / 2)), NpmALICE(m_nsubs * (10 / 2));
    vector<int> NpipALICE(m_nsubs * (10 / 2)), NpimALICE(m_nsubs * (10 / 2));
    vector<int> NKpALICE(m_nsubs * (10 / 2)), NKmALICE(m_nsubs * (10 / 2));
    vector<int> NLpALICE(m_nsubs * (10 / 2)), NLmALICE(m_nsubs * (10 / 2));

    // ALICE net-charge
    double pTmin = 0.2, pTmax = 5.0;
    vector<int> NchpALICE(m_nsubs * (10 / 2)), NchmALICE(m_nsubs * (10 / 2));

    // STAR (anti)protons
    double pTminSTAR = 0.4, pTmaxSTAR = 2.0;
    vector<int> NppSTAR(m_nsubs * (10 / 2)), NpmSTAR(m_nsubs * (10 / 2));

    // STAR off-diagonals
    double pTminSTARoff = 0.4, pTmaxSTARoff = 1.6;
    vector<int> NnetpSTARoff(m_nsubs / 2);
    vector<int> NnetkSTARoff(m_nsubs / 2);
    vector<int> NnetQSTARoff(m_nsubs / 2);
    vector<int> NnetpiSTARoff(m_nsubs / 2);

    // STAR/Kitazawa BQ
    double pTminKit = 0.4, pTmaxKit = 1.6;
    vector<int> NnetBKitaz(m_nsubs / 2);
    vector<int> NnetBKitazEta(m_nsubs / 2);
    vector<int> NnetBKitazpT(m_nsubs / 2);
    vector<int> NnetBKitazEtapT(m_nsubs / 2);
    vector<int> NnetQKitaz(m_nsubs / 2);
    vector<int> NnetQKitazEta(m_nsubs / 2);
    vector<int> NnetQKitazpT(m_nsubs / 2);
    vector<int> NnetQKitazEtapT(m_nsubs / 2);

    vector<int> NnetBKitazPrim(m_nsubs / 2);
    vector<int> NnetBKitazPrimEta(m_nsubs / 2);
    vector<int> NnetQKitazPrim(m_nsubs / 2);
    vector<int> NnetQKitazPrimEta(m_nsubs / 2);


    vector<int> NnetBKitazEtaS(m_nsubs / 2);
    vector<int> NnetQKitazEtaS(m_nsubs / 2);


    const double pTKmin = 0.2, pTKmax = 1.6;
    vector<int> NKpSTAR(m_nsubs / 2), NKmSTAR(m_nsubs / 2);
    const double pTLmin = 0.9, pTLmax = 2.0;
    vector<int> NLpSTAR(m_nsubs / 2), NLmSTAR(m_nsubs / 2);


    vector<int> NKp(m_nsubs / 2), NKm(m_nsubs / 2);
    vector<int> NLp(m_nsubs / 2), NLm(m_nsubs / 2);

    // Total net-charge
    int Nchptot = 0, Nchmtot = 0;
    int dNchdEta = 0;

    // Total protons and light nuclei
    int Nptot = 0, Ndtot = 0, Nttot = 0;

    vector<int> NppT(m_iterspT, 0), NKpT(m_iterspT, 0), NpipT(m_iterspT, 0);
    int NtotppT = 0, NtotKpT = 0, NtotpipT = 0;
    double pPt = 0., KpT = 0., pipT = 0.;

    vector<double> phisProt;

    for (const SimpleParticle& part : evt.Particles) {
      const ThermalParticle part_properties = TPS.ParticleByPDG(part.PDGID);

      // eta_s
      {
        //if (abs(part.GetEtaS()) > 5.8)
        //  std::cout << part.GetEtaS() << "   ";
        int tindeta = floor(abs(part.GetEtaS()) / m_dY);
        if (tindeta < NBps.size()) {
          if (part_properties.BaryonCharge() == 1) {
            NBps[tindeta]++;
            if (part.PDGID == 2212) {
              Npps[tindeta]++;
            }
          }

          if (part_properties.BaryonCharge() == -1) {
            NBms[tindeta]++;
            if (part.PDGID == -2212) {
              Npms[tindeta]++;
            }
          }
        }
      }

      // Y
      {
        int tindY = floor(abs(part.GetY()) / m_dY);
        if (tindY < NBpsY.size()) {
          if (part_properties.BaryonCharge() == 1) {
            NBpsY[tindY]++;
            if (part.PDGID == 2212) {
              NppsY[tindY]++;
            }
          }

          if (part_properties.BaryonCharge() == -1) {
            NBmsY[tindY]++;
            if (part.PDGID == -2212) {
              NpmsY[tindY]++;
            }
          }

          
          // Light nuclei
          if (part.PDGID == 1000010020) {
            NdpsY[tindY]++;
          }
          if (part.PDGID == -1000010020) {
            NdmsY[tindY]++;
          }
          if (part.PDGID == 1000010030) {
            NtpsY[tindY]++;
          }
          if (part.PDGID == -1000010030) {
            NtmsY[tindY]++;
          }
          if (part.PDGID == 1000020030) {
            NHe3psY[tindY]++;
          }
          if (part.PDGID == -1000020030) {
            NHe3msY[tindY]++;
          }
          if (part.PDGID == 1000020040) {
            NHe4psY[tindY]++;
          }
          if (part.PDGID == -1000020040) {
            NHe4msY[tindY]++;
          }
        }
      }

      // protons and light nuclei 4p
      {
        if (part.PDGID == 2212)
          Nptot++;
        if (part.PDGID == 1000010020)
          Ndtot++;
        if (part.PDGID == 1000010030)
          Nttot++;
      }

      // ALICE identified
      {
        int tindeta = floor(abs(part.GetEta()) / (m_dY / 10.));
        double p = part.GetP();
        if (tindeta < NppALICE.size() && p >= pmin && p <= pmax) {
          if (part.PDGID == 211)
            NpipALICE[tindeta]++;
          if (part.PDGID == -211)
            NpimALICE[tindeta]++;
          if (part.PDGID == 321)
            NKpALICE[tindeta]++;
          if (part.PDGID == -321)
            NKmALICE[tindeta]++;
          if (part.PDGID == 2212)
            NppALICE[tindeta]++;
          if (part.PDGID == -2212)
            NpmALICE[tindeta]++;
        }
      }

      // ALICE charges
      {
        int tindeta = floor(abs(part.GetEta()) / (m_dY / 10.));
        double pT = part.GetPt();
        if (tindeta < NchpALICE.size() && pT >= pTmin && pT <= pTmax) {
          if (part_properties.ElectricCharge() > 0)
            NchpALICE[tindeta]++;
          if (part_properties.ElectricCharge() < 0)
            NchmALICE[tindeta]++;
        }
      }

      // STAR identified
      {
        int tindY = floor(abs(part.GetY()) / (m_dY / 10.));
        double pT = part.GetPt();
        if (tindY < NppSTAR.size() && pT >= pTminSTAR && pT <= pTmaxSTAR) {
          if (part.PDGID == 2212)
            NppSTAR[tindY]++;
          if (part.PDGID == -2212)
            NpmSTAR[tindY]++;
        }

        if (part.PDGID == 2212 && abs(part.GetY()) < 0.5 && pT >= 0.4 && pT <= 0.8) {
          phisProt.push_back(atan2(part.py, part.px));
        }

        int tindeta = floor(abs(part.GetEta()) / (m_dY / 1.));
        if (tindeta < NnetpSTARoff.size() && pT >= pTminSTARoff && pT <= pTmaxSTARoff) {
          NnetQSTARoff[tindeta] += part_properties.ElectricCharge();
          if (part.PDGID == 321)
            NnetkSTARoff[tindeta]++;
          if (part.PDGID == -321)
            NnetkSTARoff[tindeta]--;
          if (part.PDGID == 2212)
            NnetpSTARoff[tindeta]++;
          if (part.PDGID == -2212)
            NnetpSTARoff[tindeta]--;
          if (part.PDGID == 211)
            NnetpiSTARoff[tindeta]++;
          if (part.PDGID == -211)
            NnetpiSTARoff[tindeta]--;
        }

        tindY = floor(abs(part.GetY()) / (m_dY / 1.));
        if (tindY < NnetBKitaz.size()) {
          NnetBKitaz[tindY] += part_properties.BaryonCharge();
          NnetQKitaz[tindY] += part_properties.ElectricCharge();
          if (pT >= pTminKit && pT <= pTmaxKit) {
            NnetBKitazpT[tindY] += part_properties.BaryonCharge();
            NnetQKitazpT[tindY] += part_properties.ElectricCharge();
          }
        }

        if (tindeta < NnetBKitazEta.size()) {
          NnetBKitazEta[tindeta] += part_properties.BaryonCharge();
          NnetQKitazEta[tindeta] += part_properties.ElectricCharge();
          if (pT >= pTminKit && pT <= pTmaxKit) {
            NnetBKitazEtapT[tindeta] += part_properties.BaryonCharge();
            NnetQKitazEtapT[tindeta] += part_properties.ElectricCharge();
          }
        }

        int tindetas = floor(abs(part.GetEtaS()) / (m_dY / 1.));
        if (tindetas < NnetBKitazEtaS.size()) {
          NnetBKitazEtaS[tindetas] += part_properties.BaryonCharge();
          NnetQKitazEtaS[tindetas] += part_properties.ElectricCharge();
        }
      }

      if (part_properties.ElectricCharge() > 0)
        Nchptot++;
      if (part_properties.ElectricCharge() < 0)
        Nchmtot++;

      if (part_properties.ElectricCharge() != 0 && fabs(part.GetEta()) < 0.5) {
        dNchdEta++;
      }

      // pT spectra at mid-rapidity
      if (abs(part.GetY()) < 0.5) {
        int tindpT = floor(part.GetPt() / m_dpT);
        if (tindpT < NppT.size()) {
          if (part.PDGID == 2212)
            NppT[tindpT]++;

          if (part.PDGID == 321)
            NKpT[tindpT]++;

          if (part.PDGID == 211)
            NpipT[tindpT]++;
        }

        if (part.PDGID == 2212) {
          NtotppT++;
          statspMeanpT.AddObservation(part.GetPt());
        }

        if (part.PDGID == 321) {
          NtotKpT++;
          KpT += part.GetPt();
          statsKMeanpT.AddObservation(part.GetPt());
        }

        if (part.PDGID == 211) {
          NtotpipT++;
          pipT += part.GetPt();
          statspiMeanpT.AddObservation(part.GetPt());
        }
      }

      
    }

    for (const SimpleParticle& part : evt.AllParticles) {
      if (part.PDGID == 321) {
        int tindY = floor(abs(part.GetY()) / (m_dY / 1.));
        if (tindY < NKpSTAR.size() && part.GetPt() >= pTKmin && part.GetPt() <= pTKmax)
          NKpSTAR[tindY]++;
        if (tindY < NKpSTAR.size())
          NKp[tindY]++;
      }
      if (part.PDGID == -321) {
        int tindY = floor(abs(part.GetY()) / (m_dY / 1.));
        if (tindY < NKmSTAR.size() && part.GetPt() >= pTKmin && part.GetPt() <= pTKmax)
          NKmSTAR[tindY]++;
        if (tindY < NKmSTAR.size())
          NKm[tindY]++;
      }
      if (part.PDGID == 3122) {
        int tindY = floor(abs(part.GetY()) / (m_dY / 1.));
        if (tindY < NLpSTAR.size() && part.GetPt() >= pTLmin && part.GetPt() <= pTLmax)
          NLpSTAR[tindY]++;
        if (tindY < NKpSTAR.size())
          NLp[tindY]++;
      }
      if (part.PDGID == -3122) {
        int tindY = floor(abs(part.GetY()) / (m_dY / 1.));
        if (tindY < NLmSTAR.size() && part.GetPt() >= pTLmin && part.GetPt() <= pTLmax)
          NLmSTAR[tindY]++;
        if (tindY < NKpSTAR.size())
          NLm[tindY]++;
      }

      if (part.epoch == 0) {
        const ThermalParticle part_properties = TPS.ParticleByPDG(part.PDGID);
        int tindY = floor(abs(part.GetY()) / (m_dY / 1.));
        if (tindY < NnetBKitazPrim.size()) {
          NnetBKitazPrim[tindY] += part_properties.BaryonCharge();
          NnetQKitazPrim[tindY] += part_properties.ElectricCharge();
        }
        int tindeta = floor(abs(part.GetEta()) / (m_dY / 1.));
        if (tindY < NnetBKitazPrimEta.size()) {
          NnetBKitazPrimEta[tindY] += part_properties.BaryonCharge();
          NnetQKitazPrimEta[tindY] += part_properties.ElectricCharge();
        }
      }
    }

    // Compute prefix sums to get the acceptance dependence
    vector<vector<int>*> all_numbers = {
      &NBps, &NBms, &Npps, &Npms,
      &NBpsY, &NBmsY, &NppsY, &NpmsY,
      &NppALICE, &NpmALICE, &NpipALICE, &NpimALICE,
      &NKpALICE, &NKmALICE, &NLpALICE, &NLmALICE,
      &NchpALICE, &NchmALICE, &NppSTAR, &NpmSTAR,
      &NnetQSTARoff, &NnetpSTARoff, &NnetkSTARoff, &NnetpiSTARoff,
      & NdpsY, & NdmsY, & NtpsY, &NtmsY,
      & NHe3psY,& NHe3msY,& NHe4psY,& NHe4msY,
      &NnetBKitaz, &NnetBKitazEta, &NnetBKitazEtapT,&NnetBKitazpT,
      & NnetQKitaz,& NnetQKitazEta,& NnetQKitazEtapT,& NnetQKitazpT,
      & NnetBKitazEtaS, & NnetQKitazEtaS,
       & NnetBKitazPrim,& NnetBKitazPrimEta,
      & NnetQKitazPrim,& NnetQKitazPrimEta,
      & NKpSTAR, & NKmSTAR,
      & NLpSTAR, & NLmSTAR,
      & NKp,& NKm,
      & NLp,& NLm
    };
    for (auto& nums : all_numbers) {
      for (int i = 0; i < (*nums).size(); ++i) {
        if (i > 0)
          (*nums)[i] += (*nums)[i - 1];
      }
    }


    vector<pair<vector<int>*, vector<int>*>> Npairs_Rapidity =
    {
      {&NBps, &NBms},
      {&Npps, &Npms},
      {&NBpsY, &NBmsY},
      {&NppsY, &NpmsY},
      {&NdpsY, &NdmsY},
      {&NtpsY, &NtmsY},
      {&NHe3psY, &NHe3msY},
      {&NHe4psY, &NHe4msY}
    };

    for (int in = 0; in < Npairs_Rapidity.size(); ++in) {
      auto& stats = statsRapidity[in];
      for (int isub = 0; isub < (*stats)[0].size(); ++isub) {
        (*stats)[0][isub].AddObservation(Npairs_Rapidity[in].first->operator[](isub));
        (*stats)[1][isub].AddObservation(Npairs_Rapidity[in].second->operator[](isub));
        (*stats)[2][isub].AddObservation(Npairs_Rapidity[in].first->operator[](isub) - Npairs_Rapidity[in].second->operator[](isub));
      }
    }

    for (int isub = 0; isub < twostatsB.size(); ++isub) {
      twostatsB[isub].AddObservation(NBps[isub], NBms[isub]);
      twostatsp[isub].AddObservation(Npps[isub], Npms[isub]);
      twostatsnetBp[isub].AddObservation(NBps[isub] - NBms[isub], Npps[isub] - Npms[isub]);

      statsnetB[isub].AddObservation(NBps[isub] - NBms[isub]);
      twostatsnetsumB[isub].AddObservation(NBps[isub] - NBms[isub], NBps[isub] + NBms[isub]);
    }

    for (int isub = 0; isub < twostatsnetpnetkSTAR.size(); ++isub) {
      twostatsnetpnetkSTAR[isub].AddObservation(NnetpSTARoff[isub], NnetkSTARoff[isub]);
      twostatsnetQnetkSTAR[isub].AddObservation(NnetQSTARoff[isub], NnetkSTARoff[isub]);
      twostatsnetQnetpSTAR[isub].AddObservation(NnetQSTARoff[isub], NnetpSTARoff[isub]);
      twostatsnetpinetpSTAR[isub].AddObservation(NnetpiSTARoff[isub], NnetpSTARoff[isub]);
      twostatsnetpinetkSTAR[isub].AddObservation(NnetpiSTARoff[isub], NnetkSTARoff[isub]);
    }

    for (int isub = 0; isub < twostatsBQKitazawa.size(); ++isub) {
      twostatsBQKitazawa[isub].AddObservation(NnetBKitaz[isub], NnetQKitaz[isub]);
      twostatsBQKitazawaEta[isub].AddObservation(NnetBKitazEta[isub], NnetQKitazEta[isub]);
      twostatsBQKitazawaEtaMix[isub].AddObservation(NnetBKitaz[isub], NnetQKitazEta[isub]);
      twostatsBQKitazawapT[isub].AddObservation(NnetBKitazpT[isub], NnetQKitazpT[isub]);
      twostatsBQKitazawapTEta[isub].AddObservation(NnetBKitazEtapT[isub], NnetQKitazEtapT[isub]);
      twostatsBQKitazawapTEtaMix[isub].AddObservation(NnetBKitazpT[isub], NnetQKitazEtapT[isub]);
      twostatsBQKitazawaEtaS[isub].AddObservation(NnetBKitazEtaS[isub], NnetQKitazEtaS[isub]);
      twostatsBQKitazawaPrim[isub].AddObservation(NnetBKitazPrim[isub], NnetQKitazPrim[isub]);
      twostatsBQKitazawaPrimEta[isub].AddObservation(NnetBKitazPrimEta[isub], NnetQKitazPrimEta[isub]);
      twostatsBQKitazawaPrimEtaMix[isub].AddObservation(NnetBKitazPrim[isub], NnetQKitazPrimEta[isub]);
    }

    for (int isub = 0; isub < twostatsnetsumKaon.size(); ++isub) {
      twostatsnetsumKaon[isub].AddObservation(NKpSTAR[isub] - NKmSTAR[isub], NKpSTAR[isub] + NKmSTAR[isub]);
      twostatsnetsumLambda[isub].AddObservation(NLpSTAR[isub] - NLmSTAR[isub], NLpSTAR[isub] + NLmSTAR[isub]);
      twostatsnetsumKaonNoCut[isub].AddObservation(NKp[isub] - NKm[isub], NKp[isub] + NKm[isub]);
      twostatsnetsumLambdaNoCut[isub].AddObservation(NLp[isub] - NLm[isub], NLp[isub] + NLm[isub]);
    }



    vector<pair<vector<int>*, vector<int>*>> Npairs_ALICE = {
      {&NppALICE, &NpmALICE},
      {&NKpALICE, &NKmALICE},
      {&NpipALICE, &NpimALICE},
      {&NLpALICE, &NLmALICE},
      {&NchpALICE, &NchmALICE}
    };

    for (int in = 0; in < Npairs_ALICE.size(); ++in) {
      auto& stats = statsALICE[in];
      for (int isub = 0; isub < (*stats)[0].size(); ++isub) {
        (*stats)[0][isub].AddObservation(Npairs_ALICE[in].first->operator[](isub));
        (*stats)[1][isub].AddObservation(Npairs_ALICE[in].second->operator[](isub));
        (*stats)[2][isub].AddObservation(Npairs_ALICE[in].first->operator[](isub) - Npairs_ALICE[in].second->operator[](isub));
      }
    }

    vector<pair<vector<int>*, vector<int>*>> Npairs_STAR = {
      {&NppSTAR, &NpmSTAR}
    };

    for (int in = 0; in < Npairs_STAR.size(); ++in) {
      auto& stats = statsSTAR[in];
      for (int isub = 0; isub < (*stats)[0].size(); ++isub) {
        (*stats)[0][isub].AddObservation(Npairs_STAR[in].first->operator[](isub));
        (*stats)[1][isub].AddObservation(Npairs_STAR[in].second->operator[](isub));
        (*stats)[2][isub].AddObservation(Npairs_STAR[in].first->operator[](isub) - Npairs_STAR[in].second->operator[](isub));
      }

    }

    for (int isub = 0; isub < NppSTAR.size(); ++isub) {
      twostatspSTAR[isub].AddObservation(NppSTAR[isub], NpmSTAR[isub]);
    }

    statsNpNmtot.AddObservation(Nchptot, Nchmtot);



    vector<vector<int>*> NpTs = {
      &NppT, &NKpT, &NpipT
    };

    for (int in = 0; in < NpTs.size(); ++in) {
      auto& stats = statspTs[in];
      for (int ipT = 0; ipT < (*stats).size(); ++ipT) {
        (*stats)[ipT].AddObservation(NpTs[in]->operator[](ipT));
      }
    }

    statstotppT.AddObservation(NtotppT);
    statstotKpT.AddObservation(NtotKpT);
    statstotpipT.AddObservation(NtotpipT);

    statsdNchdEta.AddObservation(dNchdEta);

    statsNpfull.AddObservation(Nptot);
    statsNdfull.AddObservation(Ndtot);
    statsNtfull.AddObservation(Nttot);

    for (int i1 = 0; i1 < (int)phisProt.size() - 1; ++i1) {
      for (int i2 = i1 + 1; i2 < phisProt.size(); ++i2) {
        double dphi = phisProt[i2] - phisProt[i1];
        if (dphi > xMath::Pi())
          dphi -= 2. * xMath::Pi();
        if (dphi < -xMath::Pi())
          dphi += 2. * xMath::Pi();
        protAngular.insert(dphi);
      }
    }
    protAngular.updateEvent();
  }

};

pair<double, double> GetSigma(SampleMoments::TwoNumberStatistics& stats) {
  vector<pair<int,int>> indices = { {1,0}, {0,1}, {2,0}, {0,2}, {1,1} };
  map<pair<int, int>, int> mapTo;
  for (int i = 0; i < indices.size(); ++i)
    mapTo[indices[i]] = i;
  vector<double> means;
  for (auto& ind : indices) {
    means.push_back(stats.GetJointMoment(ind.first, ind.second));
  }
  //{
  //  stats.GetMean1(), // k10
  //  stats.GetMean2(), // k01
  //  stats.GetJointCentralMoment(2,0), // k20
  //  stats.GetJointCentralMoment(0,2), // k02
  //  stats.GetJointCentralMoment(1,1)  // k11
  //};

  double Sigma = means[mapTo[{1, 0}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]] * (means[mapTo[{0, 2}]] - means[mapTo[{0, 1}]] * means[mapTo[{0, 1}]])
    + means[mapTo[{0, 1}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]] * (means[mapTo[{2, 0}]] - means[mapTo[{1, 0}]] * means[mapTo[{1, 0}]])
    - 2. * (means[mapTo[{1, 1}]] - means[mapTo[{1, 0}]] * means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]);

  vector<double> derivs = {
    1. / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]] * (means[mapTo[{0, 2}]] - means[mapTo[{0, 1}]] * means[mapTo[{0, 1}]])
    - means[mapTo[{1, 0}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]] * (means[mapTo[{0, 2}]] - means[mapTo[{0, 1}]] * means[mapTo[{0, 1}]])
    - 2. * means[mapTo[{0, 1}]] * means[mapTo[{1, 0}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]]
    - means[mapTo[{0, 1}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]] * (means[mapTo[{2, 0}]] - means[mapTo[{1, 0}]] * means[mapTo[{1, 0}]])
    + 2. * means[mapTo[{0, 1}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]])
    - 2. * (means[mapTo[{1, 1}]] - means[mapTo[{1, 0}]] * means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]),

    1. / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]] * (means[mapTo[{2, 0}]] - means[mapTo[{1,0}]] * means[mapTo[{1, 0}]])
    - means[mapTo[{1, 0}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]] * (means[mapTo[{0, 2}]] - means[mapTo[{0, 1}]] * means[mapTo[{0, 1}]])
    - 2. * means[mapTo[{1, 0}]] * means[mapTo[{0, 1}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]]
    - means[mapTo[{0, 1}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]] * (means[mapTo[{2, 0}]] - means[mapTo[{1, 0}]] * means[mapTo[{1, 0}]])
    + 2. * means[mapTo[{1, 0}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]])
    - 2. * (means[mapTo[{1, 1}]] - means[mapTo[{1, 0}]] * means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]),

    means[mapTo[{0, 1}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]],

    means[mapTo[{1, 0}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]],

    -2. / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]])
  };

  double SigErr = 0.0;

  for (int i = 0; i < indices.size(); ++i) {
    for (int j = 0; j < indices.size(); ++j) {
      SigErr += derivs[i] * derivs[j] * stats.GetJointMomentsSampleCovariance(indices[i].first, indices[i].second, indices[j].first, indices[j].second);
    }
  }

  return { Sigma, sqrt(SigErr) };
}

pair<double, double> GetDelta(SampleMoments::TwoNumberStatistics& stats) {
  vector<pair<int, int>> indices = { {1,0}, {0,1}, {2,0}, {0,2}, {1,1} };
  map<pair<int, int>, int> mapTo;
  for (int i = 0; i < indices.size(); ++i)
    mapTo[indices[i]] = i;
  vector<double> means;
  for (auto& ind : indices) {
    means.push_back(stats.GetJointMoment(ind.first, ind.second));
  }
  //{
  //  stats.GetMean1(), // k10
  //  stats.GetMean2(), // k01
  //  stats.GetJointCentralMoment(2,0), // k20
  //  stats.GetJointCentralMoment(0,2), // k02
  //  stats.GetJointCentralMoment(1,1)  // k11
  //};

  double Delta = -(means[mapTo[{1, 0}]] / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]] * (means[mapTo[{0, 2}]] - means[mapTo[{0, 1}]] * means[mapTo[{0, 1}]])
    - means[mapTo[{0, 1}]] / (means[mapTo[{1, 0}]] + means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]] * (means[mapTo[{2, 0}]] - means[mapTo[{1, 0}]] * means[mapTo[{1, 0}]]));

  vector<double> derivs = {
    -1. / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]] * (means[mapTo[{0, 2}]] - means[mapTo[{0, 1}]] * means[mapTo[{0, 1}]])
    + means[mapTo[{1, 0}]] / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]] * (means[mapTo[{0, 2}]] - means[mapTo[{0, 1}]] * means[mapTo[{0, 1}]])
    - 2. * means[mapTo[{0, 1}]] * means[mapTo[{1, 0}]] / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]]
    - means[mapTo[{0, 1}]] / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]] * (means[mapTo[{2, 0}]] - means[mapTo[{1, 0}]] * means[mapTo[{1, 0}]]),

    1. / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]] * (means[mapTo[{2, 0}]] - means[mapTo[{1,0}]] * means[mapTo[{1, 0}]])
    - means[mapTo[{1, 0}]] / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]] * (means[mapTo[{0, 2}]] - means[mapTo[{0, 1}]] * means[mapTo[{0, 1}]])
    + 2. * means[mapTo[{1, 0}]] * means[mapTo[{0, 1}]] / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]]
    + means[mapTo[{0, 1}]] / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]] * (means[mapTo[{2, 0}]] - means[mapTo[{1, 0}]] * means[mapTo[{1, 0}]]),

    means[mapTo[{0, 1}]] / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / means[mapTo[{1, 0}]],

    -means[mapTo[{1, 0}]] / (means[mapTo[{1, 0}]] - means[mapTo[{0, 1}]]) / means[mapTo[{0, 1}]],

    0.
  };

  double DelErr = 0.0;

  for (int i = 0; i < indices.size(); ++i) {
    for (int j = 0; j < indices.size(); ++j) {
      DelErr += derivs[i] * derivs[j] * stats.GetJointMomentsSampleCovariance(indices[i].first, indices[i].second, indices[j].first, indices[j].second);
    }
  }

  return { Delta, sqrt(DelErr) };
}

pair<double, double> GetNuDyn(SampleMoments::TwoNumberStatistics& stats) {
  vector<pair<int, int>> indices = { {1,0}, {0,1}, {2,0}, {0,2}, {1,1} };
  map<pair<int, int>, int> mapTo;
  for (int i = 0; i < indices.size(); ++i)
    mapTo[indices[i]] = i;
  vector<double> means;
  for (auto& ind : indices) {
    means.push_back(stats.GetJointMoment(ind.first, ind.second));
  }

  double nudyn = means[mapTo[{2, 0}]] / means[mapTo[{1, 0}]] / means[mapTo[{1, 0}]]
    + means[mapTo[{0, 2}]] / means[mapTo[{0, 1}]] / means[mapTo[{0, 1}]]
    - 1. / means[mapTo[{1, 0}]]
    - 1. / means[mapTo[{0, 1}]]
    - 2. * means[mapTo[{1, 1}]] / means[mapTo[{1, 0}]] / means[mapTo[{0, 1}]];

  vector<double> derivs = {
    -2. * means[mapTo[{2, 0}]] / means[mapTo[{1, 0}]] / means[mapTo[{1, 0}]] / means[mapTo[{1, 0}]] 
    + 1. / means[mapTo[{1, 0}]] / means[mapTo[{1, 0}]] 
    + 2. * means[mapTo[{1, 1}]] / means[mapTo[{1, 0}]] / means[mapTo[{1, 0}]] / means[mapTo[{0, 1}]],

    -2. * means[mapTo[{0, 2}]] / means[mapTo[{0, 1}]] / means[mapTo[{0, 1}]] / means[mapTo[{0, 1}]]
    + 1. / means[mapTo[{0, 1}]] / means[mapTo[{0, 1}]]
    + 2. * means[mapTo[{1, 1}]] / means[mapTo[{1, 0}]] / means[mapTo[{0, 1}]] / means[mapTo[{0, 1}]],

    1. / means[mapTo[{1, 0}]] / means[mapTo[{1, 0}]],

    1. / means[mapTo[{0, 1}]] / means[mapTo[{0, 1}]],

    -2. / means[mapTo[{1, 0}]] / means[mapTo[{0, 1}]]
  };

  double nudynErr = 0.0;

  for (int i = 0; i < indices.size(); ++i) {
    for (int j = 0; j < indices.size(); ++j) {
      nudynErr += derivs[i] * derivs[j] * stats.GetJointMomentsSampleCovariance(indices[i].first, indices[i].second, indices[j].first, indices[j].second);
    }
  }

  return { nudyn, sqrt(nudynErr) };
}

#endif