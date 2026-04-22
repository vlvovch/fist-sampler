#ifndef FILE_WRITER_FLUCTUATIONS_H
#define FILE_WRITER_FLUCTUATIONS_H

// Gather various statistics from events for fluctuation analysis

#include "sample-moments/include/NumberStatistics.h"
#include "sample-moments/include/TwoNumberStatistics.h"
#include "EventsProcessorFluctuations.h"

using namespace std;
using namespace thermalfist;
using namespace SampleMoments;

void WriteToFile(const string& prefix, EventsProcessorFluctuations& nstats) {
  ofstream fout;

  // Misc data
  {
    fout.open(prefix + "." + "misc.dat");
    fout << setw(15) << "nevents" << " ";
    fout << setw(15) << nstats.nevents << " ";
    fout << endl;
    fout.close();
  }

  // pT and dNdY spectra
  {
    fout.open(prefix + "." + "pT.spectra.MUSIC.dat");

    fout << setw(15) << "pT[GeV/c]" << " ";
    fout << setw(15) << "dNprot/dpT" << " ";
    fout << setw(15) << "error" << " ";
    fout << setw(15) << "dNkaon/dpT" << " ";
    fout << setw(15) << "error" << " ";
    fout << setw(15) << "dNpion/dpT" << " ";
    fout << setw(15) << "error" << " ";
    fout << endl;


    for (int ipT = 0; ipT < nstats.m_iterspT; ++ipT) {
      double pT = (0.5 + ipT) * nstats.m_dpT;
      fout << setw(15) << pT << " ";
      fout << setw(15) << nstats.statsppT[ipT].GetMean() / nstats.m_dpT << " ";
      fout << setw(15) << nstats.statsppT[ipT].GetMeanError() / nstats.m_dpT << " ";
      fout << setw(15) << nstats.statsKpT[ipT].GetMean() / nstats.m_dpT << " ";
      fout << setw(15) << nstats.statsKpT[ipT].GetMeanError() / nstats.m_dpT << " ";
      fout << setw(15) << nstats.statspipT[ipT].GetMean() / nstats.m_dpT << " ";
      fout << setw(15) << nstats.statspipT[ipT].GetMeanError() / nstats.m_dpT << " ";
      fout << endl;

    }

    fout << setw(15) << "pT-integrated" << " ";
    fout << setw(15) << nstats.statstotppT.GetMean() << " ";
    fout << setw(15) << nstats.statstotppT.GetMeanError() << " ";
    fout << setw(15) << nstats.statstotKpT.GetMean() << " ";
    fout << setw(15) << nstats.statstotKpT.GetMeanError() << " ";
    fout << setw(15) << nstats.statstotpipT.GetMean() << " ";
    fout << setw(15) << nstats.statstotpipT.GetMeanError() << " ";
    fout << endl;

    fout << setw(15) << "Mean-pT[GeV]" << " ";
    fout << setw(15) << nstats.statspMeanpT.GetMean() << " ";
    fout << setw(15) << nstats.statspMeanpT.GetMeanError() << " ";
    fout << setw(15) << nstats.statsKMeanpT.GetMean() << " ";
    fout << setw(15) << nstats.statsKMeanpT.GetMeanError() << " ";
    fout << setw(15) << nstats.statspiMeanpT.GetMean() << " ";
    fout << setw(15) << nstats.statspiMeanpT.GetMeanError() << " ";
    fout << endl;

    fout.close();


    fout.open(prefix + "." + "dNdY.MUSIC.dat");

    fout << setw(15) << "Y" << " ";
    fout << setw(15) << "dNp/dY" << " ";
    fout << setw(15) << "error" << " ";
    fout << setw(15) << "dNap/dY" << " ";
    fout << setw(15) << "error" << " ";
    fout << setw(15) << "dNetp/dY" << " ";
    fout << setw(15) << "error" << " ";
    fout << endl;

    double Nptotprev = 0., Naptotprev = 0., Nettotprev = 0.;
    double Nptotpreverr = 0., Naptotpreverr = 0., Nettotpreverr = 0.;

    for (int iY = 0; iY < nstats.statspY[0].size(); ++iY) {
      fout << setw(15) << 2. * nstats.m_dY * (1 + iY) << " ";
      double Nptot = nstats.statspY[0][iY].GetMean(), Nptoterr = nstats.statspY[0][iY].GetMeanError();
      double Naptot = nstats.statspY[1][iY].GetMean(), Naptoterr = nstats.statspY[1][iY].GetMeanError();
      double Nettot = nstats.statspY[2][iY].GetMean(), Nettoterr = nstats.statspY[2][iY].GetMeanError();
      fout << setw(15) << (Nptot - Nptotprev) / (2. * nstats.m_dY) << " ";
      fout << setw(15) << sqrt(Nptoterr * Nptoterr + Nptotpreverr * Nptotpreverr) / (2. * nstats.m_dY) << " ";
      fout << setw(15) << (Naptot - Naptotprev) / (2. * nstats.m_dY) << " ";
      fout << setw(15) << sqrt(Naptoterr * Naptoterr + Naptotpreverr * Naptotpreverr) / (2. * nstats.m_dY) << " ";
      fout << setw(15) << (Nettot - Nettotprev) / (2. * nstats.m_dY) << " ";
      fout << setw(15) << sqrt(Nettoterr * Nettoterr + Nettotpreverr * Nettotpreverr) / (2. * nstats.m_dY) << " ";
      fout << endl;

      Nptotprev = Nptot; Nptotpreverr = Nptoterr;
      Naptotprev = Naptot; Naptotpreverr = Naptoterr;
      Nettotprev = Nettot; Nettotpreverr = Nettoterr;
    }

    fout.close();
  }

  // pT integrated cumulants of baryons and protons
  {
    vector<string> names = { "baryons.dat","abaryons.dat","net-baryons.dat" };
    for (int iout = 0; iout < 3; ++iout) {
      fout.open(prefix + "." + names[iout]);
      fout << setw(15) << "dY" << " ";
      fout << setw(15) << "C1B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C2B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C2B/C1B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C3B/C1B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C4B/C1B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C1p" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C2p" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C2p/C1p" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C3p/C1p" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C4p/C1p" << " ";
      fout << setw(15) << "error" << " ";
      fout << endl;

      for (int iY = 0; iY < nstats.statsBY[iout].size(); ++iY) {
        fout << setw(15) << 2. * nstats.m_dY * (1 + iY) << " ";
        fout << setw(15) << nstats.statsBY[iout][iY].GetCumulant(1) << " ";
        fout << setw(15) << nstats.statsBY[iout][iY].GetCumulantError(1) << " ";
        fout << setw(15) << nstats.statsBY[iout][iY].GetCumulant(2) << " ";
        fout << setw(15) << nstats.statsBY[iout][iY].GetCumulantError(2) << " ";
        fout << setw(15) << nstats.statsBY[iout][iY].GetCumulantRatio(2, 1) << " ";
        fout << setw(15) << nstats.statsBY[iout][iY].GetCumulantRatioError(2, 1) << " ";
        fout << setw(15) << nstats.statsBY[iout][iY].GetCumulantRatio(3, 1) << " ";
        fout << setw(15) << nstats.statsBY[iout][iY].GetCumulantRatioError(3, 1) << " ";
        fout << setw(15) << nstats.statsBY[iout][iY].GetCumulantRatio(4, 1) << " ";
        fout << setw(15) << nstats.statsBY[iout][iY].GetCumulantRatioError(4, 1) << " ";
        fout << setw(15) << nstats.statspY[iout][iY].GetCumulant(1) << " ";
        fout << setw(15) << nstats.statspY[iout][iY].GetCumulantError(1) << " ";
        fout << setw(15) << nstats.statspY[iout][iY].GetCumulant(2) << " ";
        fout << setw(15) << nstats.statspY[iout][iY].GetCumulantError(2) << " ";
        fout << setw(15) << nstats.statspY[iout][iY].GetCumulantRatio(2, 1) << " ";
        fout << setw(15) << nstats.statspY[iout][iY].GetCumulantRatioError(2, 1) << " ";
        fout << setw(15) << nstats.statspY[iout][iY].GetCumulantRatio(3, 1) << " ";
        fout << setw(15) << nstats.statspY[iout][iY].GetCumulantRatioError(3, 1) << " ";
        fout << setw(15) << nstats.statspY[iout][iY].GetCumulantRatio(4, 1) << " ";
        fout << setw(15) << nstats.statspY[iout][iY].GetCumulantRatioError(4, 1) << " ";
        fout << endl;
      }
      fout.close();
    }
  }

  // STAR cumulants
  {
    // net-proton number cumulants
    vector<string> names = { "prot.0.4.2.0.dat","aprot.0.4.2.0.dat","net-p.0.4.2.0.dat" };
    for (int iout = 0; iout < 3; ++iout) {
      fout.open(prefix + "." + names[iout]);
      fout << setw(15) << "dY" << " ";
      fout << setw(15) << "C1" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C2" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C2/C1" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C3/C1" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C4/C1" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C3/C2" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C4/C2" << " ";
      fout << setw(15) << "error" << " ";
      fout << endl;

      for (int iY = 0; iY < nstats.statspSTAR[iout].size(); ++iY) {
        fout << setw(15) << 2. * nstats.m_dY * (1 + iY) / 10. << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulant(1) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantError(1) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulant(2) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantError(2) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantRatio(2, 1) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantRatioError(2, 1) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantRatio(3, 1) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantRatioError(3, 1) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantRatio(4, 1) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantRatioError(4, 1) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantRatio(3, 2) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantRatioError(3, 2) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantRatio(4, 2) << " ";
        fout << setw(15) << nstats.statspSTAR[iout][iY].GetCumulantRatioError(4, 2) << " ";
        fout << endl;
      }
      fout.close();
    }


    // Strongly intensive quantities for p/pbar pairs
    {
      fout.open(prefix + "." + "prot.SIQ.dat");
      fout << setw(15) << "dY" << " ";
      fout << setw(15) << "C1pp" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C1pm" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C2pp" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C2pm" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "C11" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "Sigma" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "Delta" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "nudyn" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "nudyn_crs" << " ";
      fout << setw(15) << "error" << " ";
      fout << endl;

      for (int iY = 0; iY < nstats.twostatspSTAR.size(); ++iY) {
        fout << setw(15) << 2. * nstats.m_dY * (1 + iY) / 10. << " ";
        fout << setw(15) << 2. * nstats.twostatspSTAR[iY].GetMean1() << " ";
        fout << setw(15) << 2. * nstats.twostatspSTAR[iY].GetMean1Error() << " ";
        fout << setw(15) << 2. * nstats.twostatspSTAR[iY].GetMean2() << " ";
        fout << setw(15) << 2. * nstats.twostatspSTAR[iY].GetMean2Error() << " ";
        fout << setw(15) << 2. * nstats.twostatspSTAR[iY].GetJointCentralMoment({ 2,0 }) << " ";
        fout << setw(15) << 2. * nstats.twostatspSTAR[iY].GetJointCentralMomentError({ 2,0 }) << " ";
        fout << setw(15) << 2. * nstats.twostatspSTAR[iY].GetJointCentralMoment({ 0,2 }) << " ";
        fout << setw(15) << 2. * nstats.twostatspSTAR[iY].GetJointCentralMomentError({ 0,2 }) << " ";
        fout << setw(15) << 2. * nstats.twostatspSTAR[iY].GetJointCentralMoment({ 1,1 }) << " ";
        fout << setw(15) << 2. * nstats.twostatspSTAR[iY].GetJointCentralMomentError({ 1,1 }) << " ";

        auto Sigma = GetSigma(nstats.twostatspSTAR[iY]);
        fout << setw(15) << Sigma.first << " ";
        fout << setw(15) << Sigma.second << " ";

        auto Delta = GetDelta(nstats.twostatspSTAR[iY]);
        fout << setw(15) << Delta.first << " ";
        fout << setw(15) << Delta.second << " ";

        auto nudyn = GetNuDyn(nstats.twostatspSTAR[iY]);
        fout << setw(15) << nudyn.first << " ";
        fout << setw(15) << nudyn.second << " ";

        double mult = (nstats.twostatspSTAR[iY].GetMean1() + nstats.twostatspSTAR[iY].GetMean2()) / nstats.twostatspSTAR[iY].GetMean1() / nstats.twostatspSTAR[iY].GetMean2();
        fout << setw(15) << mult * (Sigma.first - 1.) << " ";
        fout << setw(15) << mult * Sigma.second << " ";

        fout << endl;
      }
      fout.close();
    }

    // off-diagonal pkQpi as fct of deta
    // measurements/acceptance: https://arxiv.org/abs/1903.05370
    {
      fout.open(prefix + ".offdiag.pkQpi.dat");
      fout << setw(15) << "deta" << " ";
      fout << setw(15) << "pp" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "kk" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "QQ" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "pk" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "Qk" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "Qp" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "pk/k2" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "Qk/k2" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "Qp/p2" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "pipi" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "pik" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "pip" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "pik/k2" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "pip/p2" << " ";
      fout << setw(15) << "error" << " ";
      fout << endl;

      for (int iY = 0; iY < nstats.twostatsnetpnetkSTAR.size(); ++iY) {
        fout << setw(15) << 2. * nstats.m_dY * (1 + iY) << " ";
        fout << setw(15) << nstats.twostatsnetpnetkSTAR[iY].GetJointCumulant({2,0}) << " ";
        fout << setw(15) << nstats.twostatsnetpnetkSTAR[iY].GetJointCumulantError({2,0}) << " ";
        fout << setw(15) << nstats.twostatsnetpnetkSTAR[iY].GetJointCumulant({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsnetpnetkSTAR[iY].GetJointCumulantError({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsnetQnetkSTAR[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetQnetkSTAR[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetpnetkSTAR[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetpnetkSTAR[iY].GetJointCumulantError({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetQnetkSTAR[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetQnetkSTAR[iY].GetJointCumulantError({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetQnetpSTAR[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetQnetpSTAR[iY].GetJointCumulantError({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetpnetkSTAR[iY].GetJointCumulantRatio({ 1,1 }, {0,2}) << " ";
        fout << setw(15) << nstats.twostatsnetpnetkSTAR[iY].GetJointCumulantRatioError({ 1,1 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsnetQnetkSTAR[iY].GetJointCumulantRatio({ 1,1 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsnetQnetkSTAR[iY].GetJointCumulantRatioError({ 1,1 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsnetQnetpSTAR[iY].GetJointCumulantRatio({ 1,1 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsnetQnetpSTAR[iY].GetJointCumulantRatioError({ 1,1 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsnetpinetkSTAR[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetpinetkSTAR[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetpinetkSTAR[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetpinetkSTAR[iY].GetJointCumulantError({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetpinetpSTAR[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetpinetpSTAR[iY].GetJointCumulantError({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetpinetkSTAR[iY].GetJointCumulantRatio({ 1,1 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsnetpinetkSTAR[iY].GetJointCumulantRatioError({ 1,1 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsnetpinetpSTAR[iY].GetJointCumulantRatio({ 1,1 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsnetpinetpSTAR[iY].GetJointCumulantRatioError({ 1,1 }, { 0,2 }) << " ";
        fout << endl;
      }
      fout.close();
    }

    // BQ, a la https://arxiv.org/abs/2205.10030
    {
      fout.open(prefix + ".BQ.Kitazawa.dat");
      fout << setw(15) << "dY/deta/detaS" << " ";
      fout << setw(15) << "c2B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2B/c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "CBQ" << " ";
      fout << setw(15) << "error" << " ";

      fout << setw(15) << "c2B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2B/c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "CBQ" << " ";
      fout << setw(15) << "error" << " ";

      fout << setw(15) << "c2B/c2Q_mix" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "CBQ_mix" << " ";
      fout << setw(15) << "error" << " ";

      fout << setw(15) << "c2B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2B/c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "CBQ" << " ";
      fout << setw(15) << "error" << " ";

      fout << endl;

      for (int iY = 0; iY < nstats.twostatsBQKitazawa.size(); ++iY) {
        fout << setw(15) << 2. * nstats.m_dY * (1 + iY) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawa[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawa[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawa[iY].GetJointCumulant({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawa[iY].GetJointCumulantError({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawa[iY].GetJointCumulantRatio({2,0}, {0,2}) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawa[iY].GetJointCumulantRatioError({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawa[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawa[iY].GetJointCumulantError({ 1,1 }) << " ";

        fout << setw(15) << nstats.twostatsBQKitazawaEta[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEta[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEta[iY].GetJointCumulant({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEta[iY].GetJointCumulantError({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEta[iY].GetJointCumulantRatio({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEta[iY].GetJointCumulantRatioError({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEta[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEta[iY].GetJointCumulantError({ 1,1 }) << " ";

        fout << setw(15) << nstats.twostatsBQKitazawaEtaMix[iY].GetJointCumulantRatio({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEtaMix[iY].GetJointCumulantRatioError({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEtaMix[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEtaMix[iY].GetJointCumulantError({ 1,1 }) << " ";

        fout << setw(15) << nstats.twostatsBQKitazawaEtaS[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEtaS[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEtaS[iY].GetJointCumulant({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEtaS[iY].GetJointCumulantError({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEtaS[iY].GetJointCumulantRatio({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEtaS[iY].GetJointCumulantRatioError({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEtaS[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaEtaS[iY].GetJointCumulantError({ 1,1 }) << " ";
        fout << endl;
      }
      fout.close();

      fout.open(prefix + ".BQ.Kitazawa.primordial.dat");
      fout << setw(15) << "dY/deta" << " ";
      fout << setw(15) << "c2B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2B/c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "CBQ" << " ";
      fout << setw(15) << "error" << " ";

      fout << setw(15) << "c2B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2B/c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "CBQ" << " ";
      fout << setw(15) << "error" << " ";

      fout << setw(15) << "c2B/c2Q_mix" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "CBQ_mix" << " ";
      fout << setw(15) << "error" << " ";

      fout << endl;

      for (int iY = 0; iY < nstats.twostatsBQKitazawaPrim.size(); ++iY) {
        fout << setw(15) << 2. * nstats.m_dY * (1 + iY) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrim[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrim[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrim[iY].GetJointCumulant({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrim[iY].GetJointCumulantError({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrim[iY].GetJointCumulantRatio({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrim[iY].GetJointCumulantRatioError({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrim[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrim[iY].GetJointCumulantError({ 1,1 }) << " ";

        fout << setw(15) << nstats.twostatsBQKitazawaPrimEta[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrimEta[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrimEta[iY].GetJointCumulant({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrimEta[iY].GetJointCumulantError({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrimEta[iY].GetJointCumulantRatio({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrimEta[iY].GetJointCumulantRatioError({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrimEta[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrimEta[iY].GetJointCumulantError({ 1,1 }) << " ";

        fout << setw(15) << nstats.twostatsBQKitazawaPrimEtaMix[iY].GetJointCumulantRatio({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrimEtaMix[iY].GetJointCumulantRatioError({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrimEtaMix[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawaPrimEtaMix[iY].GetJointCumulantError({ 1,1 }) << " ";
        fout << endl;
      }
      fout.close();


      fout.open(prefix + ".BQ.Kitazawa.pT.dat");
      fout << setw(15) << "dY(deta)" << " ";
      fout << setw(15) << "c2B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2B/c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "CBQ" << " ";
      fout << setw(15) << "error" << " ";

      fout << setw(15) << "c2B" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2B/c2Q" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "CBQ" << " ";
      fout << setw(15) << "error" << " ";

      fout << setw(15) << "c2B/c2Q_mix" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "CBQ_mix" << " ";
      fout << setw(15) << "error" << " ";
      fout << endl;

      for (int iY = 0; iY < nstats.twostatsBQKitazawapT.size(); ++iY) {
        fout << setw(15) << 2. * nstats.m_dY * (1 + iY) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapT[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapT[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapT[iY].GetJointCumulant({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapT[iY].GetJointCumulantError({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapT[iY].GetJointCumulantRatio({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapT[iY].GetJointCumulantRatioError({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapT[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapT[iY].GetJointCumulantError({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapTEta[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapTEta[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapTEta[iY].GetJointCumulant({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapTEta[iY].GetJointCumulantError({ 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapTEta[iY].GetJointCumulantRatio({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapTEta[iY].GetJointCumulantRatioError({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapTEta[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapTEta[iY].GetJointCumulantError({ 1,1 }) << " ";

        fout << setw(15) << nstats.twostatsBQKitazawapTEtaMix[iY].GetJointCumulantRatio({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapTEtaMix[iY].GetJointCumulantRatioError({ 2,0 }, { 0,2 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapTEtaMix[iY].GetJointCumulant({ 1,1 }) << " ";
        fout << setw(15) << nstats.twostatsBQKitazawapTEtaMix[iY].GetJointCumulantError({ 1,1 }) << " ";
        fout << endl;
      }
      fout.close();
    }

    // net-kaon: https://arxiv.org/abs/1709.00773
    // net-Lambda: https://arxiv.org/abs/2001.06419
    {
      fout.open(prefix + ".KaonLambda.dat");

      fout << setw(15) << "dY" << " ";

      fout << setw(15) << "c1K" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2K" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c1K/c2K" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2K/Sk" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c3K/c1K" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c4K/c2K" << " ";
      fout << setw(15) << "error" << " ";

      fout << setw(15) << "c1L" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2L" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c1L/c2L" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2L/Sk" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c3L/c1L" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c4L/c2L" << " ";
      fout << setw(15) << "error" << " ";

      fout << endl;

      for (int iY = 0; iY < nstats.twostatsnetsumKaon.size(); ++iY) {
        fout << setw(15) << 2. * nstats.m_dY * (1 + iY) << " ";

        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulant({ 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulantError({ 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulantRatio({ 1,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulantRatioError({ 1,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulantRatio({ 2,0 }, { 0,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulantRatioError({ 2,0 }, { 0,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulantRatio({ 3,0 }, { 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulantRatioError({ 3,0 }, { 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulantRatio({ 4,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaon[iY].GetJointCumulantRatioError({ 4,0 }, { 2,0 }) << " ";

        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulant({ 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulantError({ 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulantRatio({ 1,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulantRatioError({ 1,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulantRatio({ 2,0 }, { 0,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulantRatioError({ 2,0 }, { 0,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulantRatio({ 3,0 }, { 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulantRatioError({ 3,0 }, { 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulantRatio({ 4,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambda[iY].GetJointCumulantRatioError({ 4,0 }, { 2,0 }) << " ";


        fout << endl;
      }
      fout.close();

      // Same but without pT cut
      fout.open(prefix + ".KaonLambdaNoPtCut.dat");

      fout << setw(15) << "dY" << " ";

      fout << setw(15) << "c1K" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2K" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c1K/c2K" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2K/Sk" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c3K/c1K" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c4K/c2K" << " ";
      fout << setw(15) << "error" << " ";

      fout << setw(15) << "c1L" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2L" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c1L/c2L" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c2L/Sk" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c3L/c1L" << " ";
      fout << setw(15) << "error" << " ";
      fout << setw(15) << "c4L/c2L" << " ";
      fout << setw(15) << "error" << " ";

      fout << endl;

      for (int iY = 0; iY < nstats.twostatsnetsumKaon.size(); ++iY) {
        fout << setw(15) << 2. * nstats.m_dY * (1 + iY) << " ";

        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulant({ 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulantError({ 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulantRatio({ 1,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulantRatioError({ 1,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulantRatio({ 2,0 }, { 0,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulantRatioError({ 2,0 }, { 0,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulantRatio({ 3,0 }, { 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulantRatioError({ 3,0 }, { 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulantRatio({ 4,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumKaonNoCut[iY].GetJointCumulantRatioError({ 4,0 }, { 2,0 }) << " ";

        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulant({ 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulantError({ 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulant({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulantError({ 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulantRatio({ 1,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulantRatioError({ 1,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulantRatio({ 2,0 }, { 0,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulantRatioError({ 2,0 }, { 0,1 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulantRatio({ 3,0 }, { 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulantRatioError({ 3,0 }, { 1,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulantRatio({ 4,0 }, { 2,0 }) << " ";
        fout << setw(15) << nstats.twostatsnetsumLambdaNoCut[iY].GetJointCumulantRatioError({ 4,0 }, { 2,0 }) << " ";


        fout << endl;
      }
      fout.close();
    }
  }
}

#endif