#include "Run.hpp"
#include "EMS.hpp"
#include "DrainTime.hpp"
#include "Pileup.hpp"
//#include "TMySQLServer.h"
//#include "TMySQLResult.h"
//#include "TMySQLRow.h"
#include "TPad.h"
#include "TF1.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TH2D.h"
//#include "TH3D.h"
#include <numeric>

#pragma once

/* define constants we need for later */
#define NANOSECOND .000000001
#define WINDOW 40000
#define bkgMov50ns8pe 0.10666
#define bkgMov50ns1000ns6pe 0.3
#define synthbkg_50_500_2 0.0
#define TAUN 877.7


// Definitions arranged alphabetically within their type categories.
/*-------------------- in Functions-bkg.cpp --------------------------*/
void beamOffBkg(Run* mcs1, Run* mcs2, EMS* ems0);
void beamOnBkg(Run* mcs1, Run* mcs2, EMS* ems0);
void bkgRun1Bkg(Run* mcs1, Run* mcs2, EMS* ems0);
void bkgRun2Bkg(Run* mcs1, Run* mcs2, EMS* ems0);
void longRunBkg(Run* run1, Run* run2, EMS* ems0);
void productionBkg(Run* run1, Run* run2, EMS* ems0);
void beamOffMonitor(Run* mcs1, Run* mcs2, EMS* ems0);

/*------------------- in Functions-expWeight.cpp  --------------------*/
measurement expWeightMon(TH1D* mon, double end);
//measurement expWeightMonVect(std::vector<input_t> &cts);
measurement expWeightMonVect(std::vector<input_t> &cts, double kappa);
measurement expWeightMonVectDT(std::vector<input_t> &cts, double kappa, double deadtime, double end);
measurement expWeightMonVectBkg(std::vector<input_t> &cts, double kappa, double deadtime, double bkg, double end);

/*------------------ in Functions-periodicTest.cpp -------------------*/
void protheroePeriodicDagStep(Run* runMCS1, Run* runMCS2, TH1D* histA, TH1D* histB);
void protheroePeriodicTest(Run* run, TH1D* hist);
void rayleighPeriodicTest(Run* run);
void rayleighPeriodicTestBkg(Run* run);

/*----------------------- in Functions-PMTstruct.cpp ----------------------*/
void bkgTTNE(Run* run, TH1D* ttneHist);
void fillDetHitsMoving(Run* run1, Run* run2, TH1D* dT);
void getTtne(Run* run, TH1D* ttneHist, int ch, double ts, double te);
void likelihoodStudy(Run* mcs1, Run* mcs2, TH2D* lNPE, TH2D* lmaxT);
int multiCoincCheck(Run* mcs1, Run* mcs2, TH2D* numPh, TH2D* coinL, TH2D* phByL);
//void savePMTHits(Run* run, TH2D* arrivalTime, TTree* phsTree);
void savePMTHits(Run* run1, Run* run2, TH1D* aTime1,TH1D* aTime2,TH1D* aTime3,TH1D* aTime4, TH2I* phsHist, int background);
void savePMTHits2D(Run* run1, Run* run2, TH2D* aTime1,TH2D* aTime2,TH2D* aTime3,TH2D* aTime4, bool background);
void singlePMTHits(Run* run1, Run* run2, TH1D* dT1, TH1D* dT2);
void savePMTHitsTiming(Run* mcs1, Run* mcs2, TH2D* hist2D, TH2D* asym2D, int background);
void singlePMTHitsMoving(Run* run1, Run* run2, TH1D* dT1, TH1D* dT2);
void writePMTBalance(Run* mcs1, Run* mcs2);

/*--------------------- in Functions-sumHists.cpp --------------------*/
void acSummer(Run* mcs1, Run* mcs2, TH1D* histSum);
void acSummerCoinc(Run* mcs1, Run* mcs2, TH1D* histSum, double holdT);
void peak1SummerS(Run* mcs1, Run* mcs2, TH1D* dagSum, TH1D* acSum);
void peak1SummerL(Run* mcs1, Run* mcs2, TH1D* dagSum, TH1D* acSum);
void totalSummer(Run* mcs1, Run* mcs2, TH1D* pmt1Sum, TH1D* pmt2Sum, TH1D* cSum, double holdT);
void writeFastEvts(Run* mcs1, Run* mcs2, TH1D* profile);
void writeFastBkgs(Run* mcs1, Run* mcs2, TH1D* profile);
//void acSummerAntiCoinc(Run* mcs1, Run* mcs2, TH1D* histSum, double holdT);
void bkgSummer(Run* mcs1, Run* mcs2, TH1D* pmt1Sum,TH1D* pmt2Sum);
void dipSummerL(Run* mcs1, Run* mcs2, TH1D* histSum);
void dipSummerS(Run* mcs1, Run* mcs2, TH1D* histSum, double holdT);
void monSummer(Run* mcs1, Run* mcs2, double fillEnd, int mon, TH1D* histSum);
/*void fillSummer(Run* mcs1, Run* mcs2, std::vector<double> spHits);
void fillSummerBa(Run* mcs1, double fillEnd, TH1D* histSum);
void fillSummerD1(Run* mcs1, double fillEnd, TH1D* histSum);
void fillSummerD2(Run* mcs1, double fillEnd, TH1D* histSum);
void fillSummerSp(Run* mcs1, double fillEnd, TH1D* histSum);
void fillSummerOl(Run* mcs1, double fillEnd, TH1D* histSum);*/
void longHoldSummer(Run* mcs1, Run* mcs2, TH1D* pmt1Sum, TH1D* pmt2Sum, TH1D* coincSum);

/*------------------------- in Functions.cpp -------------------------*/
int tagbitTagger(Run* run1, Run* run2, bool verbose);
double getFillEnd(Run* mcs1, Run* mcs2);
//double getDeadTimeCoinc(std::vector<coinc_t> &cts, double start, double end);
//double getDeadTimeSing(std::vector<input_t> &cts, double start, double end);
double getACUp(Run* mcs1);
std::vector<double> dagDips(Run* mcs1, Run* mcs2);
std::vector<double> dagStops(Run* mcs1, Run* mcs2, std::vector<double> dips);
std::vector<double> getLXOn(Run* mcs1, Run* mcs2, EMS* ems0);
std::vector<coinc_t> getSelfCoincs(std::vector<input_t> &cts, double initialTime, double integrateTime, int numPH);
std::vector<input_t> imposeDeadtime(std::vector<input_t> &cts, double deadtime);
std::vector<double> hMinGxHits(Run* mcs1, Run* mcs2, double fillEnd);
std::vector<input_t> removeElectricNoise_sing(std::vector<input_t> &cts, std::vector<coinc_t> &coinc, double deadt, int PE);
std::vector<coinc_t> removeElectricNoise_coinc(std::vector<coinc_t> &coinc, double deadt, int PE);
std::vector<input_t> removeHighUCN_sing(std::vector<input_t> &cts, std::vector<coinc_t> &coinc, double window, int maxPE);
std::vector<coinc_t> removeHighUCN_coinc(std::vector<coinc_t> &coinc, double window, int maxPE);
void fitBkgUnload(Run* mcs1, Run* mcs2, EMS* ems);


void effMeas(Run* mcs1, Run* mcs2);
void estimateN0(Run* run);
void extractObservables(Run* mcs1, Run* mcs2, const char* filename);
void fillingScanMeas(Run* mcs1);
//void fitFill(Run* mcs1, Run* mcs2, TH1D* fillFit);
void fitFill(Run* mcs1, Run* mcs2, const char* filename, int chan);
void fitFillFree(Run* mcs1, Run* mcs2, int chan);
void fitFillMulti(Run* mcs1, Run* mcs2);
void fitLastDip(Run* mcs1, Run* mcs2, TH1D* dipFit);
void normNByDip(Run* run1, Run* run2, EMS* ems0, const char* filename);
void normNByDipDet(Run* run1, Run* run2, EMS* ems0, const char* filename);
void normNByDipSing(Run* run1, Run* run2, EMS* ems0,const char* filename);
void normNByDipCoinc(Run* run1, Run* run2, EMS* ems0,const char* filename);
void holdingPressure(Run* run1, Run* run2, EMS* ems0);
void writeCoincHist(Run* run);
void activeCleaner(Run* run1, Run* run2, EMS* ems0, const char* filename);
int  coincidenceCompare(Run* run1, Run* run2);

/*------------------------- deprecated code --------------------------*/
//measurement normalizationLinearComb(Run* run);
//measurement integrateSP(Run* run);

//void allMonitorsAndDagger(Run* run);
//void baseAnalysis(Run* run);
//void bkgAtEnd(Run* run);
//void bkgAtEndByVdown(Run* run);
//void bkgBeforeDip(Run* run);
//void bkgByReverseDip(Run* run);
//void bkgDipTail(Run* run);
//void ctsByDip(Run* run);
//void findBkgEndRuns(Run* run);
//void getDiffFitHist(TH1* hist, TF1* fit, double& neg, double& pos);
//void getDiffAvgHist(TH1* hist, double avg, double& neg, double& pos);
//void getFillingStats(Run* run);
//void makeSummaryPlots(Run* run);
//void normNByDipByWindow(Run* run);
//void phs(Run* run);
//void saveAllPMTHits(std::vector<double> &hits);
//void sigToBkg(Run* run);
//void summarizeTTNE(Run* run);
//void summarizeDagDips(Run* run);
//void totalCoincCounts(Run* run);
