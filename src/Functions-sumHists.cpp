#include "../inc/Functions.hpp"
#include "../inc/ExpFill.hpp"
#include "TCanvas.h"
#include "TGraph.h"

//----------------------------------------------------------------------
// Unload Summers
//----------------------------------------------------------------------
// Active Cleaner summing -- finds coincidences during hold
void acSummer(Run* mcs1, Run* mcs2, TH1D* histSum) {
	
	double fillEnd = getFillEnd(mcs1, mcs2);
	// First check that we got UCN (1) during fill and (2) no beam during hold
	std::vector<input_t> foilCtsFill = mcs2->getCounts(
		[fillEnd](input_t x)->input_t{return x;},
		[fillEnd](input_t x)->bool{return x.ch == 13 && x.realtime < fillEnd;}
	);
	if(foilCtsFill.size() < 5000 || foilCtsFill.size() > 60000) {
		return;
	}
	std::vector<input_t> foilCtsHold = mcs2->getCounts(
		[fillEnd](input_t x)->input_t{return x;},
		[fillEnd](input_t x)->bool{return x.ch == 13 && x.realtime > fillEnd + 50.0;}
	);
	if(foilCtsHold.size() > 400) {
		return;
	}
	
	std::vector<input_t> acCts = mcs2->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 11;}
	);
	// Requires the "getSelfCoincs" function
	std::vector<coinc_t> coincs = getSelfCoincs(acCts, 50, 5000, 6);
			
	printf("Counts: ");
	std::for_each(coincs.begin(), coincs.end(), [](coinc_t x){
		if(x.realtime > 210 && x.realtime < 1710) {
			printf("%f,", x.realtime);
		}
	});
	printf("\n");
};

void peak1SummerS(Run* mcs1, Run* mcs2, TH1D* dagSum, TH1D* acSum) {
	
	// Load tagbits and timings
	int runNo = mcs1->getRunNo();
	// Load end of fill
	double fillEnd = getFillEnd(mcs1, mcs2);
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
	
	if(dagSteps.empty() || (dagSteps.end() <= dagSteps.begin())) { 
		fprintf(stderr, "Skipping run %d for no Dagger Movement\n", runNo);
		return;
	}
	
	// Begin hold after cleaning
	double holdStart = fillEnd + 50.0;
	double holdEnd = *dagSteps.begin(); // Beginning of counting time
	double bkgEnd = *(dagSteps.end()-1); // End of counting time
	

	// I've written two peak1Summer Histogram functions -- for short and long.
	if (holdEnd - holdStart > 21.0) {
		printf("Skipping hold, long storage time! (%d s) \n", (int)(holdEnd - holdStart));
		return;
	}
	
	if (bkgEnd - holdEnd < 40.0) {
		printf("Skipping hold, not enough time in counting period!\n");
		return;
	}
	
	Run* dagMCS; // Figure out which MCS is the dagger
	if ((mcs1->getCoincMode() == 1) || (mcs1->getCoincMode() == 2)) {
		dagMCS = mcs1;
	} else if ((mcs1->getCoincMode() == 3) || (mcs1->getCoincMode() == 4)) {
		dagMCS = mcs2;
	}
	int cAc = 1 + mcs2->getMCSOff(); // AC is mcs2 channel 1 (so offset somehow)
	
	std::vector<input_t> acCts = mcs2->getCounts(
		[holdEnd](input_t x)->input_t{input_t y = x; y.realtime -= (holdEnd - 10.0); return y;},
		[cAc, holdEnd, bkgEnd](input_t x)->bool{return x.ch == cAc && 
												x.realtime > holdEnd && x.realtime < bkgEnd;}
	);
		
	// Set AC coincidence settings (slightly laxer than two PMTs)
	double acCoinWindow = mcs1->getCoincWindow() * 2.0;
	double acPEWindow   = mcs1->getPeSumWindow() * 2.0;
	int acPESum = ceil(mcs1->getPeSum() / 2.0); 
	
	std::vector<coinc_t> acCoincs = getSelfCoincs(acCts,acCoinWindow,acPEWindow,acPESum);
	
	int numAC = 0.0;
	if (acCoincs.back().realtime > 0 && acCoincs.back().realtime < 500) {
		std::for_each(acCoincs.begin(),acCoincs.end(), [&acSum](coinc_t x){acSum->Fill(x.realtime);});
		numAC += (int)acCoincs.size();
	}
	
	std::vector<coinc_t> dagCts = dagMCS->getCoincCounts(
		[holdEnd](coinc_t x)->coinc_t{coinc_t y = x; y.realtime -= (holdEnd - 10.0); return y;},
		[holdEnd, bkgEnd](coinc_t x)->bool{return x.realtime > holdEnd && x.realtime < bkgEnd;}
	);
	
	int numD = 0.0;
	if (dagCts.back().realtime > 0 && dagCts.back().realtime < 500) {
		std::for_each(dagCts.begin(),dagCts.end(), [&dagSum](coinc_t x){dagSum->Fill(x.realtime);});
		numD += (int)dagCts.size();
	}
	
	printf("Added %d AC and %d dagger counts!\n",numAC,numD);
}
// Sum over the whole run (for imaging purposes
void totalSummer(Run* mcs1, Run* mcs2, TH1D* pmt1Sum, TH1D* pmt2Sum, TH1D* cSum, double holdT) {// Sum over all the dips of the histogram (with bkg sub.)
	
	// hard coding slop of 2.0s
	double slop = 2.0;
	int run = mcs1->getRunNo();
	if (run < 9600){
		holdT += 150.0+50.0;// 150s fill + 50s clean
	} else {
		holdT += 300.0+50.0;// 300s fill + 50s clean
	}
	// load tagbit timing
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	if(dagSteps.size() < 2) { 
		return; 
	}
	int nRuns = pmt1Sum->GetBinContent(0); // save Zeroth bin for an nruns check
	// might need to modify position of start/stop for our dips
	double cntStart = *(dagSteps.begin());
	double bkgStart = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
	double bkgEnd = *(dagSteps.end()-1);
	
	if (cntStart > holdT + slop || cntStart < holdT - slop) {
		printf("Skipping Run -- not specified holding time (%d s)!!\n", (int)(holdT - 350.0));
		printf("   This is a %d s hold.\n", (int)(cntStart - 350.0));
		return;
	}
	
	if (bkgEnd < bkgStart) {
		printf("Skipping due to background timing bug!\n");
		return;
	}
	
	printf("Start: %f stop: %f diff %f\n", cntStart, bkgStart, bkgStart-cntStart);
	// Determine which dagger pair we want to use
	Run* dagMCS;
	if (mcs1->getCoincMode() % 2) {
		dagMCS = mcs1;
	} else {
		dagMCS = mcs2;
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	// Check that we have an actual unload
	/*std::vector<coinc_t> dCts = dagMCS->getCoincCounts(
		[cntStart, bkgStart](coinc_t x)->coinc_t{x.realtime -= cntStart; return x;},
		[cntStart, bkgStart](coinc_t x)->bool{return x.realtime > cntStart && x.realtime < bkgStart;}
	);*/
	/*if(dCts.size() < 100 || dCts.empty()) {
		printf("Skipping run because there aren't enough coinc counts!\n");
		return;
	}*/
	
	// Load the counts from mcs1
	std::vector<input_t> d1Cts = dagMCS->getCounts(
		[](input_t x)->input_t{input_t y = x; return y;},
		[c1](input_t x)->bool{return x.ch == c1;}
	);	
	// Load the counts from mcs2
	std::vector<input_t> d2Cts = dagMCS->getCounts(
		[](input_t x)->input_t{input_t y = x; return y;},
		[c2](input_t x)->bool{return x.ch == c2;}
	);
	std::vector<coinc_t> cCts = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[](coinc_t x)->bool{return x.realtime > 0;}
	);
	
	if(cCts.back().realtime > 0) {
		pmt1Sum->SetBinContent(0,nRuns+1);
		pmt2Sum->SetBinContent(0,nRuns+1);
		cSum->SetBinContent(0,nRuns+1);
		printf("Added PMT 1 - %lu\n", d1Cts.size());
		printf("Added PMT 2 - %lu\n", d2Cts.size());
		std::for_each(d1Cts.begin(), d1Cts.end(), [&pmt1Sum](input_t x){pmt1Sum->Fill(x.realtime);});
		std::for_each(d2Cts.begin(), d2Cts.end(), [&pmt2Sum](input_t x){pmt2Sum->Fill(x.realtime);});
		std::for_each(cCts.begin(),  cCts.end(),  [&cSum](coinc_t x){cSum->Fill(x.realtime);});
	}
	
};


void peak1SummerL(Run* mcs1, Run* mcs2, TH1D* dagSum, TH1D* acSum) {
	
	// Load tagbits and timings
	int runNo = mcs1->getRunNo();
	// Load end of fill
	double fillEnd = getFillEnd(mcs1, mcs2);
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
	
	if(dagSteps.empty() || (dagSteps.end() <= dagSteps.begin())) { 
		fprintf(stderr, "Skipping run %d for no Dagger Movement\n", runNo);
		return;
	}
	
	// Begin hold after cleaning
	double holdStart = fillEnd + 50.0;
	double holdEnd = *dagSteps.begin(); // Beginning of counting time
	double bkgEnd = *(dagSteps.end()-1); // End of counting time
	
	
	// I've written two peak1Summer Histogram functions -- for short and long.
	if (holdEnd - holdStart < 1000.0) {
		printf("Skipping hold, short storage time! (%d s) \n", (int)(holdEnd - holdStart));
		return;
	}
	
	if (bkgEnd - holdEnd < 40.0) {
		printf("Skipping hold, not enough time in counting period!\n");
		return;
	}
	
	Run* dagMCS; // Figure out which MCS is the dagger
	if ((mcs1->getCoincMode() == 1) || (mcs1->getCoincMode() == 2)) {
		dagMCS = mcs1;
	} else if ((mcs1->getCoincMode() == 3) || (mcs1->getCoincMode() == 4)) {
		dagMCS = mcs2;
	}
	int cAc = 1 + mcs2->getMCSOff(); // AC is mcs2 channel 1 (so offset somehow)
	
	std::vector<input_t> acCts = mcs2->getCounts(
		[holdEnd](input_t x)->input_t{input_t y = x; y.realtime -= holdEnd; return y;},
		[cAc, holdEnd, bkgEnd](input_t x)->bool{return x.ch == cAc && 
												x.realtime > holdEnd && x.realtime < bkgEnd;}
	);
	// Set AC coincidence settings (slightly laxer than two PMTs)
	double acCoinWindow = mcs1->getCoincWindow() * 2.0;
	double acPEWindow   = mcs1->getPeSumWindow() * 2.0;
	int acPESum = ceil(mcs1->getPeSum() / 2.0); 
	
	std::vector<coinc_t> acCoincs = getSelfCoincs(acCts,acCoinWindow,acPEWindow,acPESum);
	
	int numAC = 0.0;
	if (acCoincs.back().realtime > 0 && acCoincs.back().realtime < 500) {
		std::for_each(acCoincs.begin(),acCoincs.end(), [&acSum](coinc_t x){acSum->Fill(x.realtime);});
		numAC += (int)acCoincs.size();
	}
	
	std::vector<coinc_t> dagCts = dagMCS->getCoincCounts(
		[holdEnd](coinc_t x)->coinc_t{coinc_t y = x; y.realtime -= holdEnd; return y;},
		[holdEnd, bkgEnd](coinc_t x)->bool{return x.realtime > holdEnd && x.realtime < bkgEnd;}
	);
	
	int numD = 0.0;
	if (dagCts.back().realtime > 0 && dagCts.back().realtime < 500) {
		std::for_each(dagCts.begin(),dagCts.end(), [&dagSum](coinc_t x){dagSum->Fill(x.realtime);});
		numD += (int)dagCts.size();
	}
	
	printf("Added %d AC and %d dagger counts!\n",numAC,numD);
}

	

// Sum active cleaner histograms across all runs of a given hold
// Coincidence mode 1 is coinc. with dagger, 2 is anticoinc. with dagger.
void acSummerCoinc(Run* mcs1, Run* mcs2, TH1D* histSum, double holdT) {
		
	// hard coding slop of 2.0s
	double slop = 2.0;	
	holdT += 300.0+50.0;// 300s fill + 50s clean
		
	// load tagbit timing
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	if(dagSteps.size() < 2) { 
		return; 
	}

	// Timing for dips
	double runStart = mcs2->getTagBitEvt(1<<1, 0.0, 1); // Start during first HGX pulse
	double cntStart = *(dagSteps.begin());
	double dipStart = *(dagSteps.begin() + 1);
	double bkgStart = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
	double bkgEnd = *(dagSteps.end()-1);
	
	if (cntStart > holdT + slop || cntStart < holdT - slop) {
		printf("Skipping Run -- not specified holding time (%d s)!!\n", (int)(holdT - 350.0));
		printf("   This is a %d s hold.\n", (int)(cntStart - 350.0));
		return;
	}
	
	if (bkgEnd < bkgStart) {
		printf("Skipping due to background timing bug!\n");
		return;
	}
	
	printf("Start: %f stop: %f diff %f\n", cntStart, dipStart, dipStart-cntStart);
	
	/*std::vector<input_t> foilCtsFill = mcs2->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 13 && x.realtime < 150.0;}
	);
	if(foilCtsFill.size() < 5000 || foilCtsFill.size() > 60000) {
		return;
	}
	std::vector<input_t> foilCtsHold = mcs2->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 13 && x.realtime > 200.0;}
	);
	if(foilCtsHold.size() > 400) {
		return;
	}
	*/
	//double cntStart = 0.0
	//double cntEnd = 
	
	/*
	std::vector<input_t> acCts = mcs2->getCoincCountsAC(
		[cntStart, bkgEnd](input_t x)->input_t{return x;},
		[cntStart, bkgStart](input_t x)->bool{return x.realtime > cntStart && x.realtime < bkgStart;}
	);
	*/

	std::vector<coinc_t> acCts = mcs2->getCoincCountsAC(
		[](coinc_t x)->coinc_t{return x;},
		[bkgEnd](coinc_t x)->bool{return x.realtime < bkgEnd;}
	);
	/*
	std::vector<input_t> acCts = mcs2->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 11;}
	);
	// Requires the "getSelfCoincs" function
	std::vector<input_t> coincs = getSelfCoincs(acCts, 50, 5000, 6);
		*/	
	/*
	if(acCts.size() < 100 || acCts.empty()) {
		printf("Skipping run because there aren't enough coinc counts!\n");
		return;
	}*/
	
	
	int numCts = 0;
	for (auto ii = acCts.begin(); ii < acCts.end()-1; ii++) {
		if (ii->realtime > cntStart && ii->realtime < dipStart) { 
			//printf("%f,",ii->realtime);
			numCts += 1;
		}
	}
	printf("Counts during dip 1: %d\n", numCts);

	
	// Look at TH1D* histSum. Extract the number and size of the bins
	// ROOT apparently never thought this would be a simple thing to do.
	int nBins = histSum->GetNbinsX();
	double binTime = (histSum->GetXaxis()->GetBinCenter(2) - histSum->GetXaxis()->GetBinCenter(1));
	
	if(acCts.back().realtime > 0 && acCts.back().realtime < 100000) {
		printf("Added - %lu\n", acCts.size());
		std::for_each(acCts.begin(), acCts.end(), [&histSum](coinc_t x){histSum->Fill(x.realtime);});
	}
	
};

// Sum the background counts to see if there's time dependence there too!
void bkgSummer(Run* mcs1, Run* mcs2, TH1D* pmt1Sum, TH1D* pmt2Sum) {
	
	int runNo = mcs1->getRunNo();
	// Load end of fill
	double fillEnd = getFillEnd(mcs1, mcs2);
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
	
	double holdT = mcs1->getCoincWindow(); // Dumb thing here but w/e
	// Begin hold after cleaning
	double holdStart = fillEnd + 50.0;
	double holdEnd = *dagSteps.begin();
	if ((holdEnd - holdStart > holdT + 1.0) || (holdEnd - holdStart < holdT - 1.0)) {
		printf("Skipping hold, not the correct holding time!");
		printf("   This is a %d s hold.\n", (int)(holdEnd - holdStart));
		return;
	}
	
	// Check if there is beam on during hold 
	std::vector<input_t> gvCts = mcs1->getCounts(
									[](input_t x)->input_t{return x;},
									[holdStart,holdEnd](input_t x)->bool{return (x.ch == 3 && x.realtime > holdStart +60.0 && x.realtime < holdEnd);}
								); 				
	if (gvCts.size() > 3100) { // Sorta arbitrarily choosing GV > 2 Hz
		printf("Skipping hold, rate is greater than 2 Hz!");
		return;	
	}
	
	// Determine which dagger pair we want to use
	Run* dagMCS;
	if ((mcs1->getCoincMode() == 1) || (mcs1->getCoincMode() == 2)) {
		dagMCS = mcs1;
	} else if ((mcs1->getCoincMode() == 3) || (mcs1->getCoincMode() == 4)) {
		dagMCS = mcs2;
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();

	// Get individual backgrounds
	double bkgStart = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
	double bkgEnd = *(dagSteps.end()-1);

	if (bkgEnd < bkgStart) {
		printf("Skipping hold, background timing bug!");
		return;
	}
	
	std::vector<input_t> bkgDag1 = dagMCS->getCounts(
			[bkgStart](input_t x)->input_t{input_t y = x; y.realtime -= (bkgStart); return y;},
			[c1,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c1 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	std::vector<input_t> bkgDag2 = dagMCS->getCounts(
			[bkgStart](input_t x)->input_t{input_t y = x; y.realtime -= (bkgStart); return y;},
			[c2,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c2 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	
	// Basic background subtraction (for scaling)
	//double bkgRate1 = (double)bkgDag1.size() / (bkgEnd - bkgStart);
	//double bkgRate2 = (double)bkgDag2.size() / (bkgEnd - bkgStart);
	
	// Look at TH1D*s longPMT(i). Extract the number and size of the bins
	// ROOT apparently never thought this would be a simple thing to do.
	//int nBins1 = pmt1Sum->GetNbinsX();
	//int nBins2 = pmt2Sum->GetNbinsX();
	/*double binTime1 = (pmt1Sum->GetXaxis()->GetBinCenter(2) - pmt1Sum->GetXaxis()->GetBinCenter(1));
	double binTime2 = (pmt2Sum->GetXaxis()->GetBinCenter(2) - pmt2Sum->GetXaxis()->GetBinCenter(1));
	
	// Create bkgSub histograms with the same axes
	TH1D* bkgSub1 = (TH1D*) pmt1Sum->Clone(); 
	TH1D* bkgSub2 = (TH1D*) pmt2Sum->Clone();
	for (int ii = 0; ii < nBins1; ii++) { // Background subtraction!
		bkgSub1->SetBinContent(ii, bkgRate1 * binTime1);
		bkgSub1->SetBinError(ii, sqrt(bkgRate1 * binTime1));
	}
	for (int ii = 0; ii < nBins2; ii++) {
		bkgSub2->SetBinContent(ii, bkgRate2 * binTime2);
		bkgSub2->SetBinError(ii, sqrt(bkgRate2 * binTime2));
	}*/
	
	// Clone pmtSums so we can add to them, making them empty
	//TH1D* pmt1Tmp = (TH1D*) pmt1Sum->Clone();
	//for (int ii = 0; ii < nBins1; ii++) { pmt1Tmp->SetBinContent(ii, 0); }
	//TH1D* pmt2Tmp = (TH1D*) pmt2Sum->Clone();
	//for (int ii = 0; ii < nBins2; ii++) { pmt2Tmp->SetBinContent(ii, 0); }
	
	// Load into temporary histograms
	int numC1 = 0.0;
	int numC2 = 0.0;
	if(bkgDag1.back().realtime > 0.0 && bkgDag1.back().realtime < bkgEnd-bkgStart) { // did we load stuff?
		std::for_each(bkgDag1.begin(), bkgDag1.end(), [&pmt1Sum](input_t x){pmt1Sum->Fill(x.realtime);}); // Note that x.realtime should be shifted by bkgStart
		numC1 += (int)bkgDag1.size();
	}
	if(bkgDag2.back().realtime > 0 && bkgDag2.back().realtime < bkgEnd-bkgStart) {
		std::for_each(bkgDag2.begin(), bkgDag2.end(), [&pmt2Sum](input_t x){pmt2Sum->Fill(x.realtime);});
		numC2 += (int)bkgDag2.size();
	}	
		
	// Re-scale pmtTmp errors
	//for (int ii = 0; ii < nBins1; ii++) {
	//	pmt1Tmp->SetBinError(ii, sqrt(pmt1Tmp->GetBinContent(ii))); // phUCN1));
	//}
	//for (int ii = 0; ii < nBins2; ii++) {
//		pmt2Tmp->SetBinError(ii, sqrt(pmt2Tmp->GetBinContent(ii))); // phUCN2));
//	}
	
	printf("Added %d counts to PMT1 and %d counts to PMT2\n",numC1,numC2);
	
	// Add to background subtracted counts
	//pmt1Tmp->Add(bkgSub1,-1.0);
	//pmt1Sum->Add(pmt1Tmp,1.0);
	//pmt2Tmp->Add(bkgSub2,-1.0);
	//pmt2Sum->Add(pmt2Tmp,1.0);
	
}

// Dip summer with no background subtraction
void dipSummerL(Run* mcs1, Run* mcs2, TH1D* histSum) {
	
	// load tagbit timing
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	if(dagSteps.size() < 2) { 
		return; 
	}

	// might need to modify position of start/stop for our dips
	double cntStart = *(dagSteps.begin());
	double bkgStart = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
	double bkgEnd = *(dagSteps.end()-1);
	
	if (cntStart < 1000 || bkgStart < 1500) {
		printf("Skipping Short Run!\n");
		return;
	}
	
	printf("Start: %f stop: %f diff %f\n", cntStart, bkgStart, bkgStart-cntStart);
	
	Run* dagMCS;
	if ((mcs1->getCoincMode() == 1) || (mcs1->getCoincMode() == 2)) {
		dagMCS = mcs1;
	} else if ((mcs1->getCoincMode() == 3) || (mcs1->getCoincMode() == 4)) {
		dagMCS = mcs2;
	}
	
	std::vector<coinc_t> dCts = dagMCS->getCoincCounts(
		[cntStart, bkgStart](coinc_t x)->coinc_t{x.realtime -= cntStart; return x;},
		[cntStart, bkgStart](coinc_t x)->bool{return x.realtime > cntStart && x.realtime < bkgStart;}
	);
	std::vector<coinc_t> bkgdCts = dagMCS->getCoincCounts(
		[bkgStart, bkgEnd](coinc_t x)->coinc_t{return x;},
		[bkgStart, bkgEnd](coinc_t x)->bool{return x.realtime > bkgStart + 5.0 && x.realtime < bkgEnd;}
	);
	if(dCts.size() < 100 || dCts.empty()) {
		return;
	}

	if(dCts.back().realtime > 0 && dCts.back().realtime < 100000) {
		printf("Added - %lu\n", dCts.size());
		printf("Bkg - %lu\n", bkgdCts.size());
		std::for_each(dCts.begin(), dCts.end(), [&histSum](coinc_t x){histSum->Fill(x.realtime);});
	}
};

// Sum over all the dips of the histogram (with bkg sub.)
void dipSummerS(Run* mcs1, Run* mcs2, TH1D* histSum, double holdT) {
	
	// hard coding slop of 2.0s
	double slop = 2.0;	
	//holdT += 300.0+50.0;// 300s fill + 50s clean
	holdT += 150.0+50.0;// 150s fill + 50s clean
	
	// load tagbit timing
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	if(dagSteps.size() < 2) { 
		return; 
	}

	// might need to modify position of start/stop for our dips
	double cntStart = *(dagSteps.begin());
	double bkgStart = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
	double bkgEnd = *(dagSteps.end()-1);
	
	if (cntStart > holdT + slop || cntStart < holdT - slop) {
		printf("Skipping Run -- not specified holding time (%d s)!!\n", (int)(holdT - 350.0));
		printf("   This is a %d s hold.\n", (int)(cntStart - 350.0));
		return;
	}
	
	if (bkgEnd < bkgStart) {
		printf("Skipping due to background timing bug!\n");
		printf("bkgStart = %f, bkgEnd = %f\n",bkgStart,bkgEnd);
		return;
	}
	// Determine which dagger pair we want to use
	Run* dagMCS;
	if ((mcs1->getCoincMode() == 1) || (mcs1->getCoincMode() == 2)) {
		dagMCS = mcs1;
	} else if ((mcs1->getCoincMode() == 3) || (mcs1->getCoincMode() == 4)) {
		dagMCS = mcs2;
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	printf("Start: %f stop: %f diff %f\n", cntStart, bkgStart, bkgStart-cntStart);
	std::vector<coinc_t> dCts = dagMCS->getCoincCounts(
		[cntStart, bkgStart](coinc_t x)->coinc_t{x.realtime -= cntStart; return x;},
		[cntStart, bkgStart](coinc_t x)->bool{return x.realtime > cntStart && x.realtime < bkgStart;}
	);
	std::vector<coinc_t> bkgdCts = dagMCS->getCoincCounts(
		[bkgStart, bkgEnd](coinc_t x)->coinc_t{return x;},
		[bkgStart, bkgEnd](coinc_t x)->bool{return x.realtime > bkgStart + 5.0 && x.realtime < bkgEnd;}
	);
	if(dCts.size() < 100 || dCts.empty()) {
		printf("Skipping run because there aren't enough coinc counts!\n");
		return;
	}

	// Basic background subtraction
	double bkgRate = (double)bkgdCts.size() / (bkgEnd - bkgStart);

	// Look at TH1D* histSum. Extract the number and size of the bins
	// ROOT apparently never thought this would be a simple thing to do.
	int nBins = histSum->GetNbinsX();
	double binTime = (histSum->GetXaxis()->GetBinCenter(2) - histSum->GetXaxis()->GetBinCenter(1));
	// Create a bkgSub histogram with the same axes
	TH1D* bkgSub = (TH1D*) histSum->Clone();  
	for (int ii = 0; ii < nBins; ii++) {
		bkgSub->SetBinContent(ii, bkgRate * binTime);
	}
	
	if(dCts.back().realtime > 0 && dCts.back().realtime < 100000) {
		printf("Added - %lu\n", dCts.size());
		printf("Bkg - %lu\n", bkgdCts.size());
		std::for_each(dCts.begin(), dCts.end(), [&histSum](coinc_t x){histSum->Fill(x.realtime);});
	}
	double binVal = 0.0;
	for (int ii = 0; ii < nBins; ii++) {
		binVal = histSum->GetBinContent(ii);
		histSum->SetBinContent(ii,binVal - bkgSub->GetBinContent(ii));	
		//histSum->Add(&bkgSub,-1.0);
	}	
};

// Monitor summing histogram
void monSummer(Run* mcs1, Run* mcs2, double countEnd, int mon, TH1D* histSum) {
	// Double checking...
	if (countEnd < 0) {return;}
	// Load Timings:
	double fillEnd = getFillEnd(mcs1, mcs2);
	if (countEnd > fillEnd + 50.) { // Check -- do we have "hold" in here?
		std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
	
		// Begin hold after cleaning, want to NOT count beam-on runs
		double holdStart = fillEnd + 50.0;
		double holdEnd = *dagSteps.begin();
		
		// For this I'm going to just force a 20s hold;
		if (holdEnd - holdStart > 25.) { return; }
		
		// Check if there is beam on during hold (really only for long but w/e)
		std::vector<input_t> gvCts = mcs1->getCounts(
							[](input_t x)->input_t{return x;},
							[holdStart,holdEnd](input_t x)->bool{return (x.ch == 3 && x.realtime >= holdStart + 60. && x.realtime < holdEnd);}
					); 
		if (gvCts.size() > 3100) { // Sorta arbitrarily choosing GV > 2 Hz
			printf("Skipping hold, rate is greater than 2 Hz!");
			return;	
		}
	}
	
	// Choose channel we want to pull from
	Run* MCS;
	// I've zero-indexed this to correspond to channels.
	if (mon < 6) { MCS = mcs1; }	
	else         { MCS = mcs2; mon -= 5; }
	
	int monReq = mon + MCS->getMCSOff();
	
	// Load counts on this channel
	std::vector<input_t> monCts = MCS->getCounts(
		[](input_t x)->input_t{input_t y = x; return y;},
		[countEnd, monReq](input_t x)->bool{return x.ch == monReq && x.realtime < countEnd;}
	);
	
	// And save!
	printf("Added %lu counts to monitor %d!\n",monCts.size(),mon);
	if (monCts.size() > 0) {
		std::for_each(monCts.begin(), monCts.end(), [&histSum](input_t x){histSum->Fill(x.realtime);});
	
		if(monCts.back().realtime > 0) {
			int nRuns = histSum->GetBinContent(0);
			histSum->SetBinContent(0,nRuns+1);
		}
	}
}

void longHoldSummer(Run* mcs1, Run* mcs2, TH1D* pmt1Sum, TH1D* pmt2Sum, TH1D* pmtCSum) {
	
	int runNo = mcs1->getRunNo();
	// Load end of fill
	double fillEnd = getFillEnd(mcs1, mcs2);
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
		
	int nRuns = pmt1Sum->GetBinContent(0); // save Zeroth bin for an nruns check
	// Begin hold after cleaning
	double holdStart = fillEnd + 50.0; // hard-coding in
	double holdEnd;
	if (dagSteps.size() > 0) { // Did we get movement?
		holdEnd = *dagSteps.begin();
	} else {
		return;
	}
	if ((holdEnd - holdStart > 1551.0) || (holdEnd - holdStart < 1549.0)) { // Force 1550s holds
		printf("Skipping run, not a 1550s hold! (holdStart = %f, holdEnd = %f)\n",holdStart,holdEnd);
		return;
	}
	
	// Check if there is beam on during hold (from 20s unload pk2)
	std::vector<input_t> gvCts = mcs1->getCounts(
									[](input_t x)->input_t{return x;},
									[holdStart,holdEnd](input_t x)->bool{return (x.ch == 3 && x.realtime > holdStart +60.0 && x.realtime < holdEnd);}
								); 				
	if (gvCts.size() > 3100) { // Arbitrarily choosing GV > 2 Hz as "beam n"
		printf("Skipping hold, rate is greater than 2 Hz! (Is beam on?)");
		return;	
	}
	
	// Determine which dagger pair we want to use
	// Can do low or high threshold
	int cMode = mcs1->getCoincMode();
	Run* dagMCS;
	// Here we're going to do 
	// Determine which dagger pair we want to use
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	// And figure out which channels to use
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();

	// Get individual backgrounds
	double bkgStart = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
	double bkgEnd = *(dagSteps.end()-1); // Should be ~50s. Don't worry about TD noise

	if (bkgEnd < bkgStart) {
		//printf("Skipping hold, background timing bug!");
		return;
	}
	
	// I attempted to scale by subtracting backgrounds;
	// I'm going to do this slightly differently.
	
	// Figure out our channels and type of background
	double bkg1Size = 0; // initialize at zero
	double bkg2Size = 0;
	double bkgCSize = 0;
		
	// Rate Dependent components
	double pmt1DT = 0;
	double pmt1PU = 0; // Stays as 0 except for integrated window!
	double pmt2DT = 0;
	double pmt2PU = 0; // Stays as 0 except for integrated window!
	double pmtCDT = 0;
	double pmtCPU = 0;
	
	std::vector<input_t> bkg1CtsT = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1, bkgStart, bkgEnd](input_t x)->bool{return (x.ch==c1 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	std::vector<input_t> bkg2CtsT = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2, bkgStart, bkgEnd](input_t x)->bool{return (x.ch==c2 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	std::vector<coinc_t> bkgCoincT = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[bkgStart, bkgEnd](coinc_t x)->bool{return (x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	
	if ((cMode == 5) || (cMode == 6)) { // Integrated window 
		
		// Load coincidence info from MCS box
		double coinWindow = dagMCS->getCoincWindow(); // Initial window
		double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
		int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
		
		std::vector<coinc_t> bkg1Cts = getSelfCoincs(bkg1CtsT, coinWindow, PEWindow, PESum);
		std::vector<coinc_t> bkg2Cts = getSelfCoincs(bkg2CtsT, coinWindow, PEWindow, PESum);
		
		// Deadtime Correction:
		pmt1DT = getDeadTimeCoinc(bkg1Cts, bkgStart, bkgEnd);
		pmt2DT = getDeadTimeCoinc(bkg2Cts, bkgStart, bkgEnd);
		pmtCDT = getDeadTimeCoinc(bkgCoincT, bkgStart, bkgEnd);
		
		// Pileup Correction: divided into pmt1 and pmt2 (and combined, duh)
		double ratePmt1  = (double)(bkg1CtsT.size())/(bkgEnd - bkgStart);
		double rateCoin1 = (double)(bkg1Cts.size()) /(bkgEnd - bkgStart);
		Pileup puObj1(dagMCS, ratePmt1, rateCoin1, bkgStart, bkgEnd);
		puObj1.LoadCoincStatsVec(bkg1Cts);
		pmt1PU = puObj1.CalculatePileup(bkg1Cts);
		
		double ratePmt2  = (double)(bkg2CtsT.size())/(bkgEnd - bkgStart);
		double rateCoin2 = (double)(bkg2Cts.size()) /(bkgEnd - bkgStart);
		Pileup puObj2(dagMCS, ratePmt2, rateCoin2, bkgStart, bkgEnd);
		puObj2.LoadCoincStatsVec(bkg2Cts);
		pmt2PU = puObj2.CalculatePileup(bkg2Cts);
		
		double ratePmtC  = (double)(bkg1CtsT.size() + bkg2CtsT.size())/(bkgEnd - bkgStart);
		double rateCoinC = (double)(bkgCoincT.size())/(bkgEnd - bkgStart);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, bkgStart, bkgEnd);
		double pmtCPU = puObjC.CalculatePileup(bkgCoincT);
		
		bkg1Size = (double)bkg1Cts.size() + pmt1DT + pmt1PU;
		bkg2Size = (double)bkg2Cts.size() + pmt2DT + pmt2PU;
		bkgCSize = (double)bkgCoincT.size() + pmtCDT + pmtCPU;
	} else if ((cMode == 7) || (cMode == 8)) { // Extra artificial deadtime
		std::vector<input_t> bkg1Cts = imposeDeadtime(bkg1CtsT,500*NANOSECOND);
		std::vector<input_t> bkg2Cts = imposeDeadtime(bkg2CtsT,500*NANOSECOND);
		
		// Deadtime Correction:
		pmt1DT = getDeadTimeSing(bkg1Cts, bkgStart, bkgEnd);
		pmt2DT = getDeadTimeSing(bkg2Cts, bkgStart, bkgEnd);
		pmtCDT = getDeadTimeCoinc(bkgCoincT, bkgStart, bkgEnd);
		
		// Pileup Correction:
		double ratePmtC  = (double)(bkg1CtsT.size() + bkg2CtsT.size())/(bkgEnd - bkgStart);
		double rateCoinC = (double)(bkgCoincT.size())/(bkgEnd - bkgStart);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, bkgStart, bkgEnd);
		double pmtCPU = puObjC.CalculatePileup(bkgCoincT);
		
		bkg1Size = (double)bkg1Cts.size() + pmt1DT + pmt1PU;
		bkg2Size = (double)bkg2Cts.size() + pmt2DT + pmt2PU;
		bkgCSize = (double)bkgCoincT.size() + pmtCDT + pmtCPU;
	} else if ((cMode == 9) || (cMode == 10)) {
		std::vector<input_t> bkg1Cts = removeElectricNoise_sing(bkg1CtsT,bkgCoincT,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<input_t> bkg2Cts = removeElectricNoise_sing(bkg2CtsT,bkgCoincT,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<coinc_t> bkgCCts = removeElectricNoise_coinc(bkgCoincT,24*NANOSECOND,dagMCS->getPeSum());
		
		// Deadtime Correction:
		pmt1DT = getDeadTimeSing(bkg1Cts, bkgStart, bkgEnd);
		pmt2DT = getDeadTimeSing(bkg2Cts, bkgStart, bkgEnd);
		pmtCDT = getDeadTimeCoinc(bkgCCts, bkgStart, bkgEnd);
		
		// Pileup Correction:
		double ratePmtC  = (double)(bkg1CtsT.size() + bkg2CtsT.size())/(bkgEnd - bkgStart);
		double rateCoinC = (double)(bkgCCts.size())/(bkgEnd - bkgStart);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, bkgStart, bkgEnd);
		double pmtCPU = puObjC.CalculatePileup(bkgCCts);
		
		bkg1Size = (double)bkg1Cts.size() + pmt1DT + pmt1PU;
		bkg2Size = (double)bkg2Cts.size() + pmt2DT + pmt2PU;
		bkgCSize = (double)bkgCoincT.size() + pmtCDT + pmtCPU;
	} else { // "Normal" running
		// Deadtime Correction:
		pmt1DT = getDeadTimeSing(bkg1CtsT, bkgStart, bkgEnd);
		pmt2DT = getDeadTimeSing(bkg2CtsT, bkgStart, bkgEnd);
		pmtCDT = getDeadTimeCoinc(bkgCoincT, bkgStart, bkgEnd);
		
		// Pileup Correction:
		double ratePmtC  = (double)(bkg1CtsT.size() + bkg2CtsT.size())/(bkgEnd - bkgStart);
		double rateCoinC = (double)(bkgCoincT.size())/(bkgEnd - bkgStart);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, bkgStart, bkgEnd);
		double pmtCPU = puObjC.CalculatePileup(bkgCoincT);
		
		bkg1Size = (double)bkg1CtsT.size() + pmt1DT + pmt1PU;
		bkg2Size = (double)bkg2CtsT.size() + pmt2DT + pmt2PU;
		bkgCSize = (double)bkgCoincT.size() + pmtCDT + pmtCPU;
	}
	
	// Basic background calculation -- want to avoid rate flurries
	double bkgRate1 = bkg1Size / (bkgEnd - bkgStart);
	double bkgRate2 = bkg2Size / (bkgEnd - bkgStart);
	double bkgRateC = bkgCSize / (bkgEnd - bkgStart);
	
	if (bkgRate1 > 200 || bkgRate2 > 200) {
		//printf("End of run background greater than 200 Hz!\n");
		return;
	}
	if (bkgRateC > 2) {
		//printf("End of run coincidence background greater than 200 Hz!\n");
		return;
	}
	// Calculate 5 sigma above the mean -- this is our expected maximum
	double rate1max = 5.0*sqrt(bkgRate1) + bkgRate1;
	double rate2max = 5.0*sqrt(bkgRate2) + bkgRate2;
	double rateCmax = 5.0*sqrt(bkgRateC) + bkgRateC;
	
	// Look at TH1D*s longPMT(i). Extract the number and size of the bins
	int nBins1 = pmt1Sum->GetNbinsX();
	int nBins2 = pmt2Sum->GetNbinsX();
	int nBinsC = pmtCSum->GetNbinsX();
	double binTime1 = (pmt1Sum->GetXaxis()->GetBinCenter(2) - pmt1Sum->GetXaxis()->GetBinCenter(1));
	double binTime2 = (pmt2Sum->GetXaxis()->GetBinCenter(2) - pmt2Sum->GetXaxis()->GetBinCenter(1));
	double binTimeC = (pmtCSum->GetXaxis()->GetBinCenter(2) - pmtCSum->GetXaxis()->GetBinCenter(1));
		
	// Clone pmtSums so we can add to them; make their clones empty
	// Ignore the first bin -- that's our counter bin
	TH1D* pmt1Tmp = (TH1D*) pmt1Sum->Clone();
	for (int ii = 1; ii < nBins1; ii++) { pmt1Tmp->SetBinContent(ii, 0); }
	pmt1Tmp->SetBinContent(0,1);
	TH1D* pmt2Tmp = (TH1D*) pmt2Sum->Clone();
	for (int ii = 1; ii < nBins2; ii++) { pmt2Tmp->SetBinContent(ii, 0); }
	pmt2Tmp->SetBinContent(0,1);
	TH1D* pmtCTmp = (TH1D*) pmtCSum->Clone();
	for (int ii = 1; ii < nBinsC; ii++) { pmtCTmp->SetBinContent(ii, 0); }
	pmtCTmp->SetBinContent(0,1);
	
	// Figure out our channels and type of background
	double pmt1Size = 0; // initialize at zero
	double pmt2Size = 0;
	double pmtCSize = 0;
	
	std::vector<input_t> pmt1CtsT = dagMCS->getCounts(
		[holdStart](input_t x)->input_t{input_t y = x; y.realtime -= holdStart; return y;},
		[c1, holdStart, holdEnd](input_t x)->bool{return (x.ch==c1 && x.realtime > holdStart && x.realtime < holdEnd);}
	);
	std::vector<input_t> pmt2CtsT = dagMCS->getCounts(
		[holdStart](input_t x)->input_t{input_t y = x; y.realtime -= holdStart; return y;},
		[c2, holdStart, holdEnd](input_t x)->bool{return (x.ch==c2 && x.realtime > holdStart && x.realtime < holdEnd);}
	);
	std::vector<coinc_t> coincT = dagMCS->getCoincCounts(
		[holdStart](coinc_t x)->coinc_t{coinc_t y = x; y.realtime -= holdStart; return y;},
		[holdStart, holdEnd](coinc_t x)->bool{return (x.realtime > holdStart && x.realtime < holdEnd);}
	);
	// Figure out our channels and type of background
		
	// Reset Rate Dependent components
	pmt1DT = 0;
	pmt1PU = 0; // Stays as 0 except for integrated window!
	pmt2DT = 0;
	pmt2PU = 0; // Stays as 0 except for integrated window!
	pmtCDT = 0;
	pmtCPU = 0;

	// Now we modify our counts based on cMode.
	if ((cMode == 5) || (cMode == 6)) { // Integrated window 
		// Load coincidence info from MCS box
		double coinWindow = dagMCS->getCoincWindow(); // Initial window
		double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
		int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
		
		std::vector<coinc_t> pmt1Cts = getSelfCoincs(pmt1CtsT, coinWindow, PEWindow, PESum);
		std::vector<coinc_t> pmt2Cts = getSelfCoincs(pmt2CtsT, coinWindow, PEWindow, PESum);
		
		// Don't worry about the RDE here... (I haven't thought about weighted PU corrections)
		if(pmt1Cts.back().realtime > 0 && pmt1Cts.back().realtime < 1551) { // This check is a little bit old
			std::for_each(pmt1Cts.begin(), pmt1Cts.end(), [&pmt1Tmp](coinc_t x){pmt1Tmp->Fill(x.realtime);});
			pmt1Size += pmt1Cts.size();
		}
		if(pmt2Cts.back().realtime > 0 && pmt2Cts.back().realtime < 1551) {
			std::for_each(pmt2Cts.begin(), pmt2Cts.end(), [&pmt2Tmp](coinc_t x){pmt2Tmp->Fill(x.realtime);});
			pmt2Size += pmt2CtsT.size();
		}
		if(coincT.back().realtime > 0 && coincT.back().realtime < 1551) { 
			std::for_each(coincT.begin(), coincT.end(), [&pmtCTmp](coinc_t x){pmtCTmp->Fill(x.realtime);});
			pmtCSize += coincT.size();
		}
		
	}
	// Now we modify our counts based on cMode.
	if ((cMode == 5) || (cMode == 6)) { // Integrated window 
		// Load coincidence info from MCS box
		double coinWindow = dagMCS->getCoincWindow(); // Initial window
		double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
		int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
		
		std::vector<coinc_t> pmt1Cts = getSelfCoincs(pmt1CtsT, coinWindow, PEWindow, PESum);
		std::vector<coinc_t> pmt2Cts = getSelfCoincs(pmt2CtsT, coinWindow, PEWindow, PESum);
		
		// Don't worry about the RDE here for now -- that's more complicated since I don't weight
		// the coincidences individually.
		if(pmt1Cts.back().realtime > 0 && pmt1Cts.back().realtime < 1551) { // This check is a little bit old
			std::for_each(pmt1Cts.begin(), pmt1Cts.end(), [&pmt1Tmp](coinc_t x){pmt1Tmp->Fill(x.realtime);});
			pmt1Size += (double)pmt1Cts.size();
		}
		if(pmt2Cts.back().realtime > 0 && pmt2Cts.back().realtime < 1551) {
			std::for_each(pmt2Cts.begin(), pmt2Cts.end(), [&pmt2Tmp](coinc_t x){pmt2Tmp->Fill(x.realtime);});
			pmt2Size += (double)pmt2CtsT.size();
		}
		if(coincT.back().realtime > 0 && coincT.back().realtime < 1551) { 
			std::for_each(coincT.begin(), coincT.end(), [&pmtCTmp](coinc_t x){pmtCTmp->Fill(x.realtime);});
			pmtCSize += (double)coincT.size();
		}
		
	} else if ((cMode == 7) || (cMode == 8)) { // Extra artificial deadtime
		std::vector<input_t> pmt1Cts = imposeDeadtime(pmt1CtsT,500*NANOSECOND);
		std::vector<input_t> pmt2Cts = imposeDeadtime(pmt2CtsT,500*NANOSECOND);
			
		if(pmt1Cts.back().realtime > 0 && pmt1Cts.back().realtime < 1551) { // This check is a little bit old
			std::for_each(pmt1Cts.begin(), pmt1Cts.end(), [&pmt1Tmp](input_t x){pmt1Tmp->Fill(x.realtime);});
			pmt1Size += (double)pmt1Cts.size();
		}
		if(pmt2Cts.back().realtime > 0 && pmt2Cts.back().realtime < 1551) {
			std::for_each(pmt2Cts.begin(), pmt2Cts.end(), [&pmt2Tmp](input_t x){pmt2Tmp->Fill(x.realtime);});
			pmt2Size += (double)pmt2CtsT.size();
		}
		if(coincT.back().realtime > 0 && coincT.back().realtime < 1551) { 
			std::for_each(coincT.begin(), coincT.end(), [&pmtCTmp](coinc_t x){pmtCTmp->Fill(x.realtime);});
			pmtCSize += (double)coincT.size();
		}
		
	} else if ((cMode == 9) || (cMode == 10)) {
		std::vector<input_t> pmt1Cts = removeElectricNoise_sing(pmt1CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<input_t> pmt2Cts = removeElectricNoise_sing(pmt2CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<coinc_t> pmtCCts = removeElectricNoise_coinc(coincT,24*NANOSECOND,dagMCS->getPeSum());
		
		if(pmt1Cts.back().realtime > 0 && pmt1Cts.back().realtime < 1551) { // This check is a little bit old
			std::for_each(pmt1Cts.begin(), pmt1Cts.end(), [&pmt1Tmp](input_t x){pmt1Tmp->Fill(x.realtime);});
			pmt1Size += (double)pmt1Cts.size();
		} 
		if(pmt2Cts.back().realtime > 0 && pmt2Cts.back().realtime < 1551) {
			std::for_each(pmt2Cts.begin(), pmt2Cts.end(), [&pmt2Tmp](input_t x){pmt2Tmp->Fill(x.realtime);});
			pmt2Size += (double)pmt2Cts.size();
		}
		if(pmtCCts.back().realtime > 0 && pmtCCts.back().realtime < 1551) { 
			std::for_each(pmtCCts.begin(), pmtCCts.end(), [&pmtCTmp](coinc_t x){pmtCTmp->Fill(x.realtime);});
			pmtCSize += (double)pmtCCts.size();
		}
			
	} else { // "Normal" running		
		if(pmt1CtsT.back().realtime > 0 && pmt1CtsT.back().realtime < 1551) { // This check is a little bit old
			std::for_each(pmt1CtsT.begin(), pmt1CtsT.end(), [&pmt1Tmp](input_t x){pmt1Tmp->Fill(x.realtime);});
			pmt1Size += (double)pmt1CtsT.size();
		}
		if(pmt2CtsT.back().realtime > 0 && pmt2CtsT.back().realtime < 1551) {
			std::for_each(pmt2CtsT.begin(), pmt2CtsT.end(), [&pmt2Tmp](input_t x){pmt2Tmp->Fill(x.realtime);});
			pmt2Size += (double)pmt2CtsT.size();
		}
		if(coincT.back().realtime > 0 && coincT.back().realtime < 1551) { 
			std::for_each(coincT.begin(), coincT.end(), [&pmtCTmp](coinc_t x){pmtCTmp->Fill(x.realtime);});
			pmtCSize += (double)coincT.size();
		}		
		
	}
	
	// Check that we're doing OK on the background sizes
	if (pmt1Size > 310000 || pmt2Size > 310000) {
		printf("Skipping hold, background rate (singles) is higher than 200 Hz!\n");
		// Side note: This is by PMT, so really 400Hz is my cutoff! That's high!
		return;
	}
	if (pmtCSize > 3100) {
		printf("Skipping hold, background rate (coinc) is higher than 2 Hz!\n");
		return;
	}
	
	// Now let's do a 5-sigma filter across the histogram:
	for (int ii = 1; ii < nBins1; ii++) {
		if (pmt1Tmp->GetBinContent(ii) > rate1max) { // Now we check if we're getting overflow stuff
			//printf("DEBUG: pmt1 bin %d, has %f (bkg ~%f)\n", ii, pmt1Tmp->GetBinContent(ii), rate1max);
			// Overflow content should just set to "average."
			pmt1Size -= pmt1Tmp->GetBinContent(ii);
			pmt1Size += bkgRate1;			
			// This is probably slightly sketchy but it should average 
			// out over hundreads of runs.
			pmt1Tmp->SetBinContent(ii,bkgRate1);
			
		}
		// And do the same for the remaining ones
		if (pmt2Tmp->GetBinContent(ii) > rate2max) {
			//printf("DEBUG: pmt2 bin %d, has %f (bkg ~%f)\n", ii, pmt2Tmp->GetBinContent(ii), rate2max);
			pmt2Size -= pmt2Tmp->GetBinContent(ii);
			pmt2Size += bkgRate2;
			pmt2Tmp->SetBinContent(ii,bkgRate2);
		}
		if (pmtCTmp->GetBinContent(ii) > rateCmax) {
			//printf("DEBUG: coinc bin %d, has %f (bkg ~%f)\n", ii, pmtCTmp->GetBinContent(ii), rateCmax);
			pmtCSize -= pmtCTmp->GetBinContent(ii);
			pmtCSize += bkgRateC;
			pmtCTmp->SetBinContent(ii,bkgRateC);
		}
	}

	// Add to our histogram	
	printf("Run %d: Added %f counts to PMT1, %f counts to PMT2, and %f to Coinc.\n",runNo, pmt1Size,pmt2Size,pmtCSize);
	
	pmt1Sum->Add(pmt1Tmp,1.0);
	pmt2Sum->Add(pmt2Tmp,1.0);
	pmtCSum->Add(pmtCTmp,1.0);
	
}

void writeFastEvts(Run* mcs1, Run* mcs2, TH1D* profile) {
	
	// Load tagbit data: FillEnd, H-GX, and Dagger Movement
	int runNo = mcs1->getRunNo();
	int cMode = mcs1->getCoincMode();
	double fillEnd = getFillEnd(mcs1, mcs2);
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2, fillEnd);
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
	std::vector<double> dagStop  = dagStops(mcs1,mcs2,dagSteps);
	// Make sure the run actually loaded...
	if(beamHits.empty() || (beamHits.end() <= beamHits.begin())) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", runNo);
		return;
	}
	if(dagSteps.empty() || (dagSteps.end() <= dagSteps.begin())) { 
		fprintf(stderr, "Skipping run %d for no Dagger Movement\n", runNo);
		return;
	}
	
	// Load our data's "run" structure
	Run* dagMCS;
	// Determine which dagger pair we want to use
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	// Get Timings	
	double holdT  = *dagSteps.begin() - fillEnd - 50.; // include cleaning.
	double unlEnd;   // Want to separate the last 50s of unload
	double countStart; // Dagger moves
	double countEnd; // TD moves
	double bkgStart; // During production running, the last 50s should have TD open, dagger down bkgs
	double bkgEnd = *(dagSteps.end()-1);
	if (dagSteps.size() > 2) { 
		countStart = *dagSteps.begin(); // just so I don't have to retype it
		countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
		unlEnd   = countEnd - 50;
		bkgStart = mcs1->getTagBitEvt(1<<4, *(dagSteps.end()-2), 1) + 5.0;// 1<<4 is Cat Door, add extra 5 seconds to be safe
	} else {
		fprintf(stderr, "Skipping run %d because not enough steps!\n",runNo);
	}
	
	if (holdT < 0.) {
		fprintf(stderr,"Skipping run %d for negative holding time\n",runNo);
		return;
	}
	
	if (bkgStart < 5.0 || bkgStart >= bkgEnd) { // tag = -1.0 for errors (added 5s for stability.)
		fprintf(stderr, "Skipping run %d for background timing glitch\n", runNo);
		return;
	}
	if (unlEnd < 0 || unlEnd < countStart) {
		fprintf(stderr,"Skipping run %f for unload timing glitch\n", runNo);
	}
	printf("Timings: %f,%f,%f,%f,%f\n",countStart,countEnd,unlEnd,bkgStart,bkgEnd);
	// Let's just take all electric counts and "blob" them together, we can separate later
	// Load all counts between counting start and end of background
	std::vector<coinc_t> cts = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[countStart, bkgEnd](coinc_t x)->bool{return (x.realtime > countStart && x.realtime < bkgEnd);}
	);

	std::vector<coinc_t> ctsF = removeElectricNoise_coinc(cts,24*NANOSECOND,dagMCS->getPeSum());
	printf("Counts: %lu,%lu\n",cts.size(),ctsF.size());
	// Really we're looking at "Removed" coincidences (those are the bad ones)
	std::vector<coinc_t> removed;
		
	// for this script I only care about removed coincidences.
	// Sort by time first
	unsigned long pk1Cts = 0;
	unsigned long pk2Cts = 0;
	unsigned long pk3Cts = 0;
	unsigned long unlCts = 0;
	unsigned long bkgCts = 0;
		
	if (cts.size() > ctsF.size()) {
		// Now need to figure out where "removed" coincidences occur.
		size_t fInd = 0;
		for (size_t cInd = 0; cInd < cts.size(); cInd++) { // Loop through indices
			if (fInd != ctsF.size() - 1) { // Basically a check if the last Coinc is bad
				if (cts.at(cInd).realtime == ctsF.at(fInd).realtime) { // Check timing
					fInd += 1; // Go on to next index
					continue;
				}
			}
			// We skipped one! Don't increment fInd, add to next vector.
			printf("time: %f\n",cts.at(cInd).realtime);
			removed.push_back(cts.at(cInd));
			// And add to timing counts.
			if (cts.at(cInd).realtime>=  *(dagSteps.begin()) && cts.at(cInd).realtime < *(dagSteps.begin()+1)) {
				pk1Cts += 1;
			} else if (cts.at(cInd).realtime >=  *(dagSteps.begin()+1) && cts.at(cInd).realtime < *(dagSteps.begin()+2)) {
				pk2Cts += 1;
			} else if (cts.at(cInd).realtime >=  *(dagSteps.begin()+2) && cts.at(cInd).realtime < unlEnd) {	
				pk3Cts += 1;
			} else if (cts.at(cInd).realtime >=  unlEnd && cts.at(cInd).realtime < countEnd) {	
				unlCts += 1;
			} else if (cts.at(cInd).realtime >=  bkgStart && cts.at(cInd).realtime < bkgEnd) {
				bkgCts += 1;
			}	
		}
	}
	
	// Now we can (1) add to our histogram and (2) write out our data
	int nRuns = profile->GetBinContent(0); // save zeroth bin to get nruns 
	profile->SetBinContent(0,nRuns+1);
	if (removed.size() > 0) {
		std::for_each(removed.begin(),removed.end(),[&profile,countStart](coinc_t x) {profile->Fill(x.realtime-countStart);});
	}
	FILE* outfile;
	outfile = fopen("electricFilter.csv","a");
	printf("Run %d (hold %f):\n",runNo,holdT);
	printf("   %lu,%lu,%lu,%lu,%lu\n",pk1Cts,pk2Cts,pk3Cts,unlCts,bkgCts);
	fprintf(outfile,"%d,%f,%lu,%lu,%lu,%lu,%lu\n",runNo,holdT,pk1Cts,pk2Cts,pk3Cts,unlCts,bkgCts);
	fclose(outfile);
	
}

void writeFastBkgs(Run* mcs1, Run* mcs2, TH1D* profile) {
	
	// Load tagbit data: FillEnd, H-GX, and Dagger Movement
	int runNo = mcs1->getRunNo();
	int cMode = mcs1->getCoincMode();
	
	double dagMvs[5]  = {0.0, 250.0, 500.0, 750.0, 1000.0};
	double stepTime;
	for (int i = 0; i < 3; i++) {
		stepTime = mcs2->getTagBitEvt(1<<0, dagMvs[i] + 2.0, 1);
		if (stepTime < 0 ) {
			printf("Error! One or more dagger steps does not exist! Returning.\n");
			return;
		}
		dagMvs[i+1] = stepTime;
	}
	
	// Load our data's "run" structure
	Run* dagMCS;
	// Determine which dagger pair we want to use
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	std::vector<input_t> gv = mcs1->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return (x.ch == 3 && x.realtime < 300);}
	);
	if (gv.size() > 600) {
		printf("Beam is on!\n");
		return;
	}
	
	// Let's just take all electric counts and "blob" them together, we can separate later
	// Load all counts between counting start and end of background
	std::vector<coinc_t> cts = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[](coinc_t x)->bool{return (x.realtime > 0. && x.realtime < 1000.);}
	);
	std::vector<coinc_t> ctsF = removeElectricNoise_coinc(cts,24*NANOSECOND,dagMCS->getPeSum());
	
	// Really we're looking at "Removed" coincidences (those are the bad ones)
	std::vector<coinc_t> removed;
		
	// for this script I only care about removed coincidences.
	// Sort by time first
	unsigned long cts38 = 0;
	unsigned long cts49 = 0;
	unsigned long cts25 = 0;
	unsigned long cts1 = 0;
			
	if (cts.size() > ctsF.size()) {
		// Now need to figure out where "removed" coincidences occur.
		size_t fInd = 0;
		for (size_t cInd = 0; cInd < cts.size(); cInd++) { // Loop through indices
			if (fInd != ctsF.size() - 1) { // Basically a check if the last Coinc is bad
				if (cts.at(cInd).realtime == ctsF.at(fInd).realtime) { // Check timing
					fInd += 1; // Go on to next index
					continue;
				}
			}
			// We skipped one! Don't increment fInd, add to next vector.
			//printf("time: %f\n",cts.at(cInd).realtime);
			removed.push_back(cts.at(cInd));
			// And add to timing counts.
			if (cts.at(cInd).realtime>=  dagMvs[0] && cts.at(cInd).realtime < dagMvs[1]) {
				cts38 += 1;
			} else if (cts.at(cInd).realtime >=  dagMvs[1] && cts.at(cInd).realtime < dagMvs[2]) {
				cts49 += 1;
			} else if (cts.at(cInd).realtime >=  dagMvs[2] && cts.at(cInd).realtime < dagMvs[3]) {	
				cts25 += 1;
			} else if (cts.at(cInd).realtime >=  dagMvs[3] && cts.at(cInd).realtime < dagMvs[4]) {	
				cts1 += 1;
			}	
		}
	}
	
	// Now we can (1) add to our histogram and (2) write out our data
	int nRuns = profile->GetBinContent(0); // save zeroth bin to get nruns 
	profile->SetBinContent(0,nRuns+1);
	if (removed.size() > 0) {
		std::for_each(removed.begin(),removed.end(),[&profile](coinc_t x) {profile->Fill(x.realtime);});
	}
	FILE* outfile;
	outfile = fopen("electricFilter.csv","a");
	printf("Run %d (hold %f):\n",runNo,-2.);
	printf("   %lu,%lu,%lu,%lu,%lu\n",cts38,cts49,cts25,cts1,0);
	fprintf(outfile,"%d,%f,%lu,%lu,%lu,%lu,%lu\n",runNo,cts38,cts49,cts25,cts1,0);
	fclose(outfile);
	
}
