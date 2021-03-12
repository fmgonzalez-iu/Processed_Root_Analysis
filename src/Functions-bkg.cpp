#include "../inc/Functions.hpp"
#include "../inc/ExpFill.hpp"
#include "TCanvas.h"
#include "TGraph.h"

/* define constants we need for later */
#define NANOSECOND .000000001
#define bkgMov50ns8pe 0.10666
#define bkgMov50ns1000ns6pe 0.3
#define synthbkg_50_500_2 0.0
#define TAUN 877.7

/*----------------------------------------------------------------------------
	Editor: Frank M. Gonzalez
	
	This file contains a bunch of background summing functions.
	
	Condensed here to make it easier(?) to modify.
	* 
	* 
	TODO: Make the "types" of coincidence easier to modify
------------------------------------------------------------------------------*/

/* dedicated beam-off background run. */
/* Should be 4 step with dagger at 380, 490, 250, and 10 for 250s each. */
void beamOffBkg(Run* mcs1, Run* mcs2, EMS* ems0) {
		
	int runNo = mcs1->getRunNo();
	
	// No beam on! 
	// This means we have to do a separate dagger movement finder, can't use dagDips.
	// Hardcoding in heights and dagger guesses. -- now fixed, remember to update python(!)
	int heights[4]     = {380, 490, 250, 10};
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
		
	// Output a .csv file that we will use in background_analyzer.py
	FILE* outfile;
	//outfile = fopen("bkgBeamOff.csv", "a"); // DEBUG
	outfile = fopen("bkgAll.csv", "a");
	printf("Outputting data to bkgAll.csv!\n");
	
	FILE* outfileRDE;
	outfileRDE = fopen("bkgRDE.csv","a");
	printf("Outputting RDE data to bkgRDE.csv!\n");
	// Loop through each dagger step and find when photon events hit
	for(int i = 1; i < 5; i++) {
		double startTime = dagMvs[i-1] + 10.0;
		double endTime = dagMvs[i]- 10.0;
				
		double pmt1Size = 0; // initialize at zero
		double pmt2Size = 0;
		double pmtCSize = 0;
		
		// Rate Dependent components
		double pmt1DT = 0;
		double pmt1PU = 0; // Stays as 0 except for integrated window!
		double pmt2DT = 0;
		double pmt2PU = 0; // Stays as 0 except for integrated window!
		double pmtCDT = 0;
		double pmtCPU = 0;
		
		std::vector<input_t> pmt1CtsT = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[c1, startTime, endTime](input_t x)->bool{return (x.ch==c1 && x.realtime >= startTime && x.realtime < endTime);}
		);
		std::vector<input_t> pmt2CtsT = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[c2, startTime, endTime](input_t x)->bool{return (x.ch==c2 && x.realtime >= startTime && x.realtime < endTime);}
		);
		std::vector<coinc_t> coincT = dagMCS->getCoincCounts(
			[](coinc_t x)->coinc_t{return x;},
			[startTime, endTime](coinc_t x)->bool{return (x.realtime >= startTime && x.realtime < endTime);}
		);
		
		// Now we modify our counts based on cMode.
		if ((cMode == 5) || (cMode == 6)) { // Integrated window 
			// Load coincidence info from MCS box
			double coinWindow = dagMCS->getCoincWindow(); // Initial window
			double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
			int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
			
			std::vector<coinc_t> pmt1Cts = getSelfCoincs(pmt1CtsT, coinWindow, PEWindow, PESum);
			std::vector<coinc_t> pmt2Cts = getSelfCoincs(pmt2CtsT, coinWindow, PEWindow, PESum);
			
			// Deadtime Correction:
			pmt1DT = getDeadTimeCoinc(pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeCoinc(pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction: divided into pmt1 and pmt2 (and combined, duh)
			double ratePmt1  = (double)(pmt1CtsT.size())/(endTime - startTime);
			double rateCoin1 = (double)(pmt1Cts.size()) /(endTime - startTime);
			Pileup puObj1(dagMCS, ratePmt1, rateCoin1, startTime, endTime);
			puObj1.LoadCoincStatsVec(pmt1Cts);
			pmt1PU = puObj1.CalculatePileup(pmt1Cts);
			
			double ratePmt2  = (double)(pmt2CtsT.size())/(endTime - startTime);
			double rateCoin2 = (double)(pmt2Cts.size()) /(endTime - startTime);
			Pileup puObj2(dagMCS, ratePmt2, rateCoin2, startTime, endTime);
			puObj2.LoadCoincStatsVec(pmt2Cts);
			pmt2PU = puObj2.CalculatePileup(pmt2Cts);
			
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		} else if ((cMode == 7) || (cMode == 8)) { // Extra artificial deadtime
			std::vector<input_t> pmt1Cts = imposeDeadtime(pmt1CtsT,500*NANOSECOND);
			std::vector<input_t> pmt2Cts = imposeDeadtime(pmt2CtsT,500*NANOSECOND);
			
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		} else if ((cMode == 9) || (cMode == 10)) {
			std::vector<input_t> pmt1Cts = removeElectricNoise_sing(pmt1CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
			std::vector<input_t> pmt2Cts = removeElectricNoise_sing(pmt2CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
			std::vector<coinc_t> pmtCCts = removeElectricNoise_coinc(coincT,24*NANOSECOND,dagMCS->getPeSum());
			
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(pmtCCts,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1Cts.size() + pmt2Cts.size())/(endTime - startTime);
			double rateCoinC = (double)(pmtCCts.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;		
		} else { // "Normal" running
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1CtsT,startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2CtsT,startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1CtsT.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2CtsT.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		}

		// Get PMT temp averages (it'll return -1 if no EMS data can be found)
		double temp1 = ems0->getEMSAvg(3311,0,startTime,endTime);
		double temp2 = ems0->getEMSAvg(3311,1,startTime,endTime);

		// Incorporate Pileup into p

		printf("\n BACKGROUNDS --- Run: %d\nheight: %d\nstartTime: %f\nendTime: %f\npmt1cts: %f\npmt2Cts: %f\ncoincCts: %f\nPMT Temps: %f, %f\n",
				runNo, heights[i-1], startTime, endTime, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
		fprintf(outfile,"%d,%d,%f,%f,%f,%f,%f,%f,%f\n",
				runNo, heights[i-1], startTime, endTime, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
		printf("RDE: \npmt1DT: %f\npmt1PU: %f\npmt2DT: %f\npmt2PU: %f\npmtCDT: %f\npmtCPU: %f\n",
				runNo, heights[i-1],pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
		fprintf(outfileRDE,"%d,%d,%f,%f,%f,%f,%f,%f\n",
				runNo, heights[i-1],pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
	}
	fclose(outfile);
	fclose(outfileRDE);
	return;
}

/* Should be 4 step with dagger at 380, 490, 250, and 10 for 250s each. */
void beamOffMonitor(Run* mcs1, Run* mcs2, EMS* ems0) {	
	int runNo = mcs1->getRunNo();
	
	// No beam on! 
	// This means we have to do a separate dagger movement finder, can't use dagDips.
	// Hardcoding in heights and dagger guesses.
	int heights[4]     = {380, 490, 250, 10};
	double dagMvs[5]  = {0.0, 250.0, 500.0, 750.0, 1000.0};
	double stepTime;
	// Just to make sure the run actually exists
	for (int i = 0; i < 3; i++) {
		stepTime = mcs2->getTagBitEvt(1<<0, dagMvs[i] + 2.0, 1);
		if (stepTime < 0 ) {
			printf("Error! One or more dagger steps does not exist! Returning.\n");
			return;
		}
		dagMvs[i+1] = stepTime;
	}
		
	// Output a .csv file that we will use in background_analyzer.py
	// I'm setting this as the same format of the previous detectorAndBkg.csv
	FILE* outfileDet;
	outfileDet = fopen("detectorAndBkg_BeamOff.csv","a");
	printf("Outputting monitor+background data to detectorAndBkg_BeamOff.csv!\n");
	
	// Loop through each dagger step and find when photon events hit
	for(long i = 1; i < 5; i++) {
		double startTime = dagMvs[i-1] + 10.0; // timings
		double endTime = dagMvs[i]- 10.0;
		std::vector<long> mon;  //Bakground of monitor counts
		
		for (int jj = 1; jj < 11; jj++) { // loop across the monitors
			int detCh;
			Run* dagMCS;
			// need to convert detectors 1-10 to the right MCS box
			if (jj < 6) {
				dagMCS = mcs1;
				detCh = jj + dagMCS->getMCSOff();
			} else {
				dagMCS = mcs2;
				detCh = jj-5 + dagMCS->getMCSOff();
			}
			std::vector<input_t> detCts = dagMCS->getCounts(
				[](input_t x)->input_t{return x;},
				[detCh, startTime,endTime](input_t x)->bool{return x.ch == detCh && 
															(x.realtime >= startTime && x.realtime < endTime);}
			);
			printf("%d, %lu\n",detCh,detCts.size());
			mon.push_back(detCts.size()); // Raw monitor counts
		}
		
		// I'm setting this as the same format of the previous detectorAndBkg.csv
		// I've converted tHold --> dip and bkgTime --> tElapse
		printf("\nMonitor Data -- Run: %d\nDip: %d\ntElapse: %f\n", 
				runNo, i, endTime - startTime);
		fprintf(outfileDet,"%d,%d,%f",runNo, i, endTime - startTime);
		// Assume the uncertainty is 0 (since it's just counting)
		for (size_t jt = 0; jt < mon.size(); jt++ ) {
			printf("Monitor %d: %ld (%u)\n",jt,mon[jt],0);
			fprintf(outfileDet,",%ld,%u",mon[jt],0);
		}
		fprintf(outfileDet,"\n");
	}
	fclose(outfileDet);
	return;
}

/* Normal holds but with GV closed */
void beamOnBkg(Run* mcs1, Run* mcs2, EMS* ems0) {
	
	// Load tagbits: TD and H-GX
	int runNo = mcs1->getRunNo();
	
	// Figure out when the dagger moves
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
	
	// If there are no dagger movements, skip the run (it was probably killed)
	if (dagSteps.end() == dagSteps.begin()) {
		fprintf(stderr, "Skipping run %d for no dagger movement!\n", runNo);
		return;
	} 
	
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
	
	// Output a .csv file that we can use
	FILE* outfile;
	//outfile = fopen("bkgByDip.csv", "a");
	outfile = fopen("bkgAll.csv", "a");
	printf("Outputting data to bkgByDip.csv\n!");	
	
	FILE* outfileRDE;
	outfileRDE = fopen("bkgRDE.csv","a");
	printf("Outputting RDE data to bkgRDE.csv!\n");
	
	int i = 0;
	int heights[3] = {380,250,10};

	// Load coincidence info for PMT (if doing single PMT integrated window)
	double coinWindow = dagMCS->getCoincWindow(); // Initial window
	double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
	int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
	
	// Loop through each dagger step and find when photon events hit
	for(auto stepIt = dagSteps.begin()+1; stepIt < dagSteps.end(); stepIt++) {
		double startTime = *(stepIt-1);
		double endTime = *(stepIt);
						
		// In the actual data, we cut off the last 50s for backgrounds.
		// We should do the same here
		if (endTime == *(dagSteps.end()-1)) { 
			// trap door moving time
			endTime = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
			// Some backgrounds have no trap door movement, block last 50s
			if (endTime < 0.0) { 
				endTime = *stepIt - 50.0;
				printf("Warning! Using last 50s for backgrounds! TD doesn't move!\n");
			}
		}		
		if (endTime < 0.0 || startTime >= endTime) {
			fprintf(stderr, "Skipping run %d for background timing glitch\n", runNo);
			return;
		}
		
		double pmt1Size = 0; // initialize at zero
		double pmt2Size = 0;
		double pmtCSize = 0;
		
		// Rate Dependent components
		double pmt1DT = 0;
		double pmt1PU = 0; // Stays as 0 except for integrated window!
		double pmt2DT = 0;
		double pmt2PU = 0; // Stays as 0 except for integrated window!
		double pmtCDT = 0;
		double pmtCPU = 0;
		
		std::vector<input_t> pmt1CtsT = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[c1, startTime, endTime](input_t x)->bool{return (x.ch==c1 && x.realtime > startTime && x.realtime < endTime);}
		);
		std::vector<input_t> pmt2CtsT = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[c2, startTime, endTime](input_t x)->bool{return (x.ch==c2 && x.realtime > startTime && x.realtime < endTime);}
		);
		std::vector<coinc_t> coincT = dagMCS->getCoincCounts(
			[](coinc_t x)->coinc_t{return x;},
			[startTime, endTime](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
		);
		// Now we modify our counts based on cMode.
		if ((cMode == 5) || (cMode == 6)) { // Integrated window 
			// Load coincidence info from MCS box
			double coinWindow = dagMCS->getCoincWindow(); // Initial window
			double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
			int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
			
			std::vector<coinc_t> pmt1Cts = getSelfCoincs(pmt1CtsT, coinWindow, PEWindow, PESum);
			std::vector<coinc_t> pmt2Cts = getSelfCoincs(pmt2CtsT, coinWindow, PEWindow, PESum);
			
			// Deadtime Correction:
			pmt1DT = getDeadTimeCoinc(pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeCoinc(pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction: divided into pmt1 and pmt2 (and combined, duh)
			double ratePmt1  = (double)(pmt1CtsT.size())/(endTime - startTime);
			double rateCoin1 = (double)(pmt1Cts.size()) /(endTime - startTime);
			Pileup puObj1(dagMCS, ratePmt1, rateCoin1, startTime, endTime);
			puObj1.LoadCoincStatsVec(pmt1Cts);
			pmt1PU = puObj1.CalculatePileup(pmt1Cts);
			
			double ratePmt2  = (double)(pmt2CtsT.size())/(endTime - startTime);
			double rateCoin2 = (double)(pmt2Cts.size()) /(endTime - startTime);
			Pileup puObj2(dagMCS, ratePmt2, rateCoin2, startTime, endTime);
			puObj2.LoadCoincStatsVec(pmt2Cts);
			pmt2PU = puObj2.CalculatePileup(pmt2Cts);
			
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		} else if ((cMode == 7) || (cMode == 8)) { // Extra artificial deadtime
			std::vector<input_t> pmt1Cts = imposeDeadtime(pmt1CtsT,500*NANOSECOND);
			std::vector<input_t> pmt2Cts = imposeDeadtime(pmt2CtsT,500*NANOSECOND);
			
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		} else if ((cMode == 9) || (cMode == 10)) {
			std::vector<input_t> pmt1Cts = removeElectricNoise_sing(pmt1CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
			std::vector<input_t> pmt2Cts = removeElectricNoise_sing(pmt2CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
			std::vector<coinc_t> pmtCCts = removeElectricNoise_coinc(coincT,24*NANOSECOND,dagMCS->getPeSum());
			
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(pmtCCts,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1Cts.size() + pmt2Cts.size())/(endTime - startTime);
			double rateCoinC = (double)(pmtCCts.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;		
		} else { // "Normal" running
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1CtsT,startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2CtsT,startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1CtsT.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2CtsT.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		}
		// Get PMT temp averages (it'll return -1 if no EMS data can be found)
		double temp1 = ems0->getEMSAvg(3311,0,startTime,endTime);
		double temp2 = ems0->getEMSAvg(3311,1,startTime,endTime);
			
		printf("\n BACKGROUNDS --- Run: %d\nheight: %d\nstartTime: %f\nendTime: %f\npmt1cts: %f\npmt2Cts: %f\ncoincCts: %f\nPMT Temps: %f, %f\n",
			runNo, heights[i], startTime, endTime, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
		fprintf(outfile,"%d,%d,%f,%f,%f,%f,%f,%f,%f\n",
			runNo, heights[i], startTime, endTime, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
		printf("RDE: \npmt1DT: %f\npmt1PU: %f\npmt2DT: %f\npmt2PU: %f\npmtCDT: %f\npmtCPU: %f\n",
				runNo, heights[i],pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
		fprintf(outfileRDE,"%d,%d,%f,%f,%f,%f,%f,%f\n",
				runNo, heights[i],pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
		i=i+1;
	}
	
	// Now that we've done the dips (height dependent, do the other spots
	longRunBkg(mcs1,mcs2,ems0);
	
	// "Normal" background spots (last 50s and hold)
	double bkgStart = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
	double bkgEnd = *(dagSteps.end()-1);
	
	if (bkgStart < 0.0) {
		bkgStart = bkgEnd - 50.0;
		printf("Hard coding in last 50s as normal background!\n");
	}
	
	if (bkgStart < 0.0 || bkgEnd < 0.0){
		fprintf(stderr," Warning! Backgrounds are being weird! Quitting run %d!\n",runNo);
		return;
	}

	double pmt1Size = 0; // initialize at zero
	double pmt2Size = 0;
	double pmtCSize = 0;
	
	// Rate Dependent components
	double pmt1DT = 0;
	double pmt1PU = 0; // Stays as 0 except for integrated window!
	double pmt2DT = 0;
	double pmt2PU = 0; // Stays as 0 except for integrated window!
	double pmtCDT = 0;
	double pmtCPU = 0;

	std::vector<input_t> pmt1CtsT = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1, bkgStart, bkgEnd](input_t x)->bool{return (x.ch==c1 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	std::vector<input_t> pmt2CtsT = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2, bkgStart, bkgEnd](input_t x)->bool{return (x.ch==c2 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	std::vector<coinc_t> coincT = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[bkgStart, bkgEnd](coinc_t x)->bool{return (x.realtime > bkgStart && x.realtime < bkgEnd);}
	);

	// Now we modify our counts based on cMode.
	if ((cMode == 5) || (cMode == 6)) { // Integrated window 
		// Load coincidence info from MCS box
		double coinWindow = dagMCS->getCoincWindow(); // Initial window
		double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
		int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
		
		std::vector<coinc_t> pmt1Cts = getSelfCoincs(pmt1CtsT, coinWindow, PEWindow, PESum);
		std::vector<coinc_t> pmt2Cts = getSelfCoincs(pmt2CtsT, coinWindow, PEWindow, PESum);
		
		// Deadtime Correction:
		pmt1DT = getDeadTimeCoinc(pmt1Cts, bkgStart, bkgEnd);
		pmt2DT = getDeadTimeCoinc(pmt2Cts, bkgStart, bkgEnd);
		pmtCDT = getDeadTimeCoinc(coincT,  bkgStart, bkgEnd);
		
		// Pileup Correction: divided into pmt1 and pmt2 (and combined, duh)
		double ratePmt1  = (double)(pmt1CtsT.size())/(bkgEnd - bkgStart);
		double rateCoin1 = (double)(pmt1Cts.size()) /(bkgEnd - bkgStart);
		Pileup puObj1(dagMCS, ratePmt1, rateCoin1, bkgStart, bkgEnd);
		puObj1.LoadCoincStatsVec(pmt1Cts);
		pmt1PU = puObj1.CalculatePileup(pmt1Cts);
		
		double ratePmt2  = (double)(pmt2CtsT.size())/(bkgEnd - bkgStart);
		double rateCoin2 = (double)(pmt2Cts.size()) /(bkgEnd - bkgStart);
		Pileup puObj2(dagMCS, ratePmt2, rateCoin2, bkgStart, bkgEnd);
		puObj2.LoadCoincStatsVec(pmt2Cts);
		pmt2PU = puObj2.CalculatePileup(pmt2Cts);
		
		double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(bkgEnd - bkgStart);
		double rateCoinC = (double)(coincT.size())/(bkgEnd - bkgStart);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, bkgStart, bkgEnd);
		double pmtCPU = puObjC.CalculatePileup(coincT);
		
		pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
		pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
		pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
	} else if ((cMode == 7) || (cMode == 8)) { // Extra artificial deadtime
		std::vector<input_t> pmt1Cts = imposeDeadtime(pmt1CtsT,500*NANOSECOND);
		std::vector<input_t> pmt2Cts = imposeDeadtime(pmt2CtsT,500*NANOSECOND);
		
		// Add RDE
		// Deadtime Correction:
		pmt1DT = getDeadTimeSing( pmt1Cts, bkgStart, bkgEnd);
		pmt2DT = getDeadTimeSing( pmt2Cts, bkgStart, bkgEnd);
		pmtCDT = getDeadTimeCoinc(coincT,  bkgStart, bkgEnd);
		
		// Pileup Correction
		double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(bkgEnd - bkgStart);
		double rateCoinC = (double)(coincT.size())/(bkgEnd - bkgStart);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, bkgStart, bkgEnd);
		double pmtCPU = puObjC.CalculatePileup(coincT);
		
		pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
		pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
		pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
	} else if ((cMode == 9) || (cMode == 10)) {
		std::vector<input_t> pmt1Cts = removeElectricNoise_sing(pmt1CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<input_t> pmt2Cts = removeElectricNoise_sing(pmt2CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<coinc_t> pmtCCts = removeElectricNoise_coinc(coincT,24*NANOSECOND,dagMCS->getPeSum());
		
		// Add RDE
		// Deadtime Correction:
		pmt1DT = getDeadTimeSing( pmt1Cts, bkgStart, bkgEnd);
		pmt2DT = getDeadTimeSing( pmt2Cts, bkgStart, bkgEnd);
		pmtCDT = getDeadTimeCoinc(pmtCCts, bkgStart, bkgEnd);
		
		// Pileup Correction
		double ratePmtC  = (double)(pmt1Cts.size() + pmt2Cts.size())/(bkgEnd - bkgStart);
		double rateCoinC = (double)(pmtCCts.size())/(bkgEnd - bkgStart);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, bkgStart, bkgEnd);
		double pmtCPU = puObjC.CalculatePileup(coincT);
		
		pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
		pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
		pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;		
	} else { // "Normal" running
		// Add RDE
		// Deadtime Correction:
		pmt1DT = getDeadTimeSing(pmt1CtsT, bkgStart, bkgEnd);
		pmt2DT = getDeadTimeSing(pmt2CtsT, bkgStart, bkgEnd);
		pmtCDT = getDeadTimeCoinc(coincT,  bkgStart, bkgEnd);
		
		// Pileup Correction
		double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(bkgEnd - bkgStart);
		double rateCoinC = (double)(coincT.size())/(bkgEnd - bkgStart);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, bkgStart, bkgEnd);
		double pmtCPU = puObjC.CalculatePileup(coincT);
		
		pmt1Size = ((double)pmt1CtsT.size()) + pmt1DT + pmt1PU;
		pmt2Size = ((double)pmt2CtsT.size()) + pmt2DT + pmt2PU;
		pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
	}
	
	// Get PMT temp averages (it'll return -1 if no EMS data can be found)
	double temp1 = ems0->getEMSAvg(3311,0,bkgStart,bkgEnd);
	double temp2 = ems0->getEMSAvg(3311,1,bkgStart,bkgEnd);
	
	// Output -- hardcoded in final height 10	
	printf("\n BACKGROUNDS --- Run: %d\nheight: %d\nstartTime: %f\nendTime: %f\npmt1cts: %f\npmt2Cts: %f\ncoincCts: %f\nPMT Temps: %f, %f\n",
		runNo, 10, bkgStart, bkgEnd, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
	fprintf(outfile,"%d,%d,%f,%f,%f,%f,%f,%f,%f\n",
		runNo, 10, bkgStart, bkgEnd, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
	printf("RDE: \npmt1DT: %f\npmt1PU: %f\npmt2DT: %f\npmt2PU: %f\npmtCDT: %f\npmtCPU: %f\n",
			runNo, 10, pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
	fprintf(outfileRDE,"%d,%d,%f,%f,%f,%f,%f,%f\n",
			runNo, 10, pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
			
	fclose(outfile);
	fclose(outfileRDE);
	return;
}

/* dedicated beam-on background run, 300s fill then 100s at each height */
void bkgRun1Bkg(Run* mcs1, Run* mcs2, EMS* ems0) {
	
	// Load tagbit data: FillEnd, H-GX, and Dagger Movement
	int runNo = mcs1->getRunNo();
	double fillEnd = getFillEnd(mcs1, mcs2);
	
	// Should be 2 dagger steps + end of fill.
	std::vector<double> dagMvs = dagDips(mcs1,mcs2);
	// Make sure the run actually loaded...
	//if(dagMvs.empty() || (dagMvs.end() <= dagMvs.begin())) { 
	if (dagMvs.end() == dagMvs.begin()) { 
		fprintf(stderr, "Skipping run %d for no Dagger Movement\n", runNo);
		return;
	}
	// Push fillEnd to the front of a separate dagSteps array (with guesses).
	double dagSteps[4] = {fillEnd, fillEnd+100.0, fillEnd+200.0, fillEnd+300.0};
	//for (int i = 0; i < dagMvs.size(); i++) {
	//	dagSteps[i+1] = dagMvs[i];
	//}
	
	// 300s "fill" followed by dagger at at 490000 -> 380000 -> 10000 for 100s each
	int heights[3] = {490, 380, 10};
	
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
	
	// Output a .csv file that we will use in NormAnalyzer-Coinc.py
	FILE* outfile;
	//outfile = fopen("bkgRunBeamOn.csv", "a");
	outfile = fopen("bkgAll.csv", "a");
	printf("Outputting data to bkgRunBeamOn.csv\n!");
		
	FILE* outfileRDE;
	outfileRDE = fopen("bkgRDE.csv","a");
	printf("Outputting RDE data to bkgRDE.csv!\n");
	// Loop through each dagger step and find when photon events hit
	//for(auto stepIt = dagSteps.begin()+1; stepIt < dagSteps.end(); stepIt++) {
	for(int i = 1; i < 4; i++) {
		double startTime = dagSteps[i-1];
		double endTime = dagSteps[i];
		
		// Figure out our channels and type of background
		double pmt1Size = 0; // initialize at zero
		double pmt2Size = 0;
		double pmtCSize = 0;
		
		// Rate Dependent components
		double pmt1DT = 0;
		double pmt1PU = 0; // Stays as 0 except for integrated window!
		double pmt2DT = 0;
		double pmt2PU = 0; // Stays as 0 except for integrated window!
		double pmtCDT = 0;
		double pmtCPU = 0;
		
		std::vector<input_t> pmt1CtsT = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[c1, startTime, endTime](input_t x)->bool{return (x.ch==c1 && x.realtime > startTime && x.realtime < endTime);}
		);
		std::vector<input_t> pmt2CtsT = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[c2, startTime, endTime](input_t x)->bool{return (x.ch==c2 && x.realtime > startTime && x.realtime < endTime);}
		);
		std::vector<coinc_t> coincT = dagMCS->getCoincCounts(
			[](coinc_t x)->coinc_t{return x;},
			[startTime, endTime](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
		);
		// Now we modify our counts based on cMode.
		if ((cMode == 5) || (cMode == 6)) { // Integrated window 
			// Load coincidence info from MCS box
			double coinWindow = dagMCS->getCoincWindow(); // Initial window
			double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
			int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
			
			std::vector<coinc_t> pmt1Cts = getSelfCoincs(pmt1CtsT, coinWindow, PEWindow, PESum);
			std::vector<coinc_t> pmt2Cts = getSelfCoincs(pmt2CtsT, coinWindow, PEWindow, PESum);
			
			// Deadtime Correction:
			pmt1DT = getDeadTimeCoinc(pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeCoinc(pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction: divided into pmt1 and pmt2 (and combined, duh)
			double ratePmt1  = (double)(pmt1CtsT.size())/(endTime - startTime);
			double rateCoin1 = (double)(pmt1Cts.size()) /(endTime - startTime);
			Pileup puObj1(dagMCS, ratePmt1, rateCoin1, startTime, endTime);
			puObj1.LoadCoincStatsVec(pmt1Cts);
			pmt1PU = puObj1.CalculatePileup(pmt1Cts);
			
			double ratePmt2  = (double)(pmt2CtsT.size())/(endTime - startTime);
			double rateCoin2 = (double)(pmt2Cts.size()) /(endTime - startTime);
			Pileup puObj2(dagMCS, ratePmt2, rateCoin2, startTime, endTime);
			puObj2.LoadCoincStatsVec(pmt2Cts);
			pmt2PU = puObj2.CalculatePileup(pmt2Cts);
			
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		} else if ((cMode == 7) || (cMode == 8)) { // Extra artificial deadtime
			std::vector<input_t> pmt1Cts = imposeDeadtime(pmt1CtsT,500*NANOSECOND);
			std::vector<input_t> pmt2Cts = imposeDeadtime(pmt2CtsT,500*NANOSECOND);
			
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		} else if ((cMode == 9) || (cMode == 10)) {
			std::vector<input_t> pmt1Cts = removeElectricNoise_sing(pmt1CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
			std::vector<input_t> pmt2Cts = removeElectricNoise_sing(pmt2CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
			std::vector<coinc_t> pmtCCts = removeElectricNoise_coinc(coincT,24*NANOSECOND,dagMCS->getPeSum());
			
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(pmtCCts,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1Cts.size() + pmt2Cts.size())/(endTime - startTime);
			double rateCoinC = (double)(pmtCCts.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;		
		} else { // "Normal" running
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1CtsT,startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2CtsT,startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1CtsT.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2CtsT.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		}
		// Get PMT temp averages (it'll return -1 if no EMS data can be found)
		double temp1 = ems0->getEMSAvg(3311,0,startTime,endTime);
		double temp2 = ems0->getEMSAvg(3311,1,startTime,endTime);
				
		printf("\n BACKGROUNDS --- Run: %d\nheight: %d\nstartTime: %f\nendTime: %f\npmt1cts: %f\npmt2Cts: %f\ncoincCts: %f\nPMT Temps: %f, %f\n",
				runNo, heights[i-1], startTime, endTime, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
		fprintf(outfile,"%d,%d,%f,%f,%f,%f,%f,%f,%f\n",
				runNo, heights[i-1], startTime, endTime, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
		printf("RDE: \npmt1DT: %f\npmt1PU: %f\npmt2DT: %f\npmt2PU: %f\npmtCDT: %f\npmtCPU: %f\n",
				runNo, heights[i-1],pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
		fprintf(outfileRDE,"%d,%d,%f,%f,%f,%f,%f,%f\n",
				runNo, heights[i-1],pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
	}
	fclose(outfile);
	fclose(outfileRDE);
	return;
}

/* Inverted dedicated background run */
void bkgRun2Bkg(Run* mcs1, Run* mcs2, EMS* ems0) {
	
	// Load tagbit data: FillEnd, H-GX, and Dagger Movement
	int runNo = mcs1->getRunNo();
	double fillEnd = getFillEnd(mcs1, mcs2);
	
	// Should be 2 dagger steps + end of fill.
	std::vector<double> dagMvs = dagDips(mcs1,mcs2);
	// Make sure the run actually loaded...
	//if(dagMvs.empty() || (dagMvs.end() <= dagMvs.begin())) { 
	if(dagMvs.end() == dagMvs.begin()) { 
		fprintf(stderr, "Skipping run %d for no Dagger Movement\n", runNo);
		return;
	}
	// Push fillEnd to the front of a separate dagSteps array (with guesses).
	double dagSteps[4] = {fillEnd, fillEnd+100.0, fillEnd+200.0, fillEnd+300.0};
	//for (int i = 0; i < dagMvs.size(); i++) {
	//	dagSteps[i+1] = dagMvs[i];
	//}
	
	// 300s "fill" followed by dagger at at 10000 -> 380000 -> 490000 for 100s each
	int heights[3] = {10, 380, 490};
	
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
	
	// Output a .csv file that we will use in NormAnalyzer-Coinc.py
	FILE* outfile;
	outfile = fopen("bkgAll.csv", "a");
	//outfile = fopen("bkgRunBeamOn.csv", "a");
	printf("Outputting data to bkgRunBeamOn.csv\n!");
		
	FILE* outfileRDE;
	outfileRDE = fopen("bkgRDE.csv","a");
	printf("Outputting RDE data to bkgRDE.csv!\n");
	
	// Loop through each dagger step and find when photon events hit
	//for(auto stepIt = dagSteps.begin()+1; stepIt < dagSteps.end(); stepIt++) {
	for(int i = 1; i < 4; i++) {
		double startTime = dagSteps[i-1];
		double endTime = dagSteps[i];
		
		
		// Figure out our channels and type of background
		double pmt1Size = 0; // initialize at zero
		double pmt2Size = 0;
		double pmtCSize = 0;
		
		// Rate Dependent components
		double pmt1DT = 0;
		double pmt1PU = 0; // Stays as 0 except for integrated window!
		double pmt2DT = 0;
		double pmt2PU = 0; // Stays as 0 except for integrated window!
		double pmtCDT = 0;
		double pmtCPU = 0;
		
		std::vector<input_t> pmt1CtsT = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[c1, startTime, endTime](input_t x)->bool{return (x.ch==c1 && x.realtime > startTime && x.realtime < endTime);}
		);
		std::vector<input_t> pmt2CtsT = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[c2, startTime, endTime](input_t x)->bool{return (x.ch==c2 && x.realtime > startTime && x.realtime < endTime);}
		);
		std::vector<coinc_t> coincT = dagMCS->getCoincCounts(
			[](coinc_t x)->coinc_t{return x;},
			[startTime, endTime](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
		);
		// Now we modify our counts based on cMode.
		if ((cMode == 5) || (cMode == 6)) { // Integrated window 
			// Load coincidence info from MCS box
			double coinWindow = dagMCS->getCoincWindow(); // Initial window
			double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
			int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
			
			std::vector<coinc_t> pmt1Cts = getSelfCoincs(pmt1CtsT, coinWindow, PEWindow, PESum);
			std::vector<coinc_t> pmt2Cts = getSelfCoincs(pmt2CtsT, coinWindow, PEWindow, PESum);
			
			// Deadtime Correction:
			pmt1DT = getDeadTimeCoinc(pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeCoinc(pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction: divided into pmt1 and pmt2 (and combined, duh)
			double ratePmt1  = (double)(pmt1CtsT.size())/(endTime - startTime);
			double rateCoin1 = (double)(pmt1Cts.size()) /(endTime - startTime);
			Pileup puObj1(dagMCS, ratePmt1, rateCoin1, startTime, endTime);
			puObj1.LoadCoincStatsVec(pmt1Cts);
			pmt1PU = puObj1.CalculatePileup(pmt1Cts);
			
			double ratePmt2  = (double)(pmt2CtsT.size())/(endTime - startTime);
			double rateCoin2 = (double)(pmt2Cts.size()) /(endTime - startTime);
			Pileup puObj2(dagMCS, ratePmt2, rateCoin2, startTime, endTime);
			puObj2.LoadCoincStatsVec(pmt2Cts);
			pmt2PU = puObj2.CalculatePileup(pmt2Cts);
			
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		} else if ((cMode == 7) || (cMode == 8)) { // Extra artificial deadtime
			std::vector<input_t> pmt1Cts = imposeDeadtime(pmt1CtsT,500*NANOSECOND);
			std::vector<input_t> pmt2Cts = imposeDeadtime(pmt2CtsT,500*NANOSECOND);
			
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		} else if ((cMode == 9) || (cMode == 10)) {
			std::vector<input_t> pmt1Cts = removeElectricNoise_sing(pmt1CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
			std::vector<input_t> pmt2Cts = removeElectricNoise_sing(pmt2CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
			std::vector<coinc_t> pmtCCts = removeElectricNoise_coinc(coincT,24*NANOSECOND,dagMCS->getPeSum());
			
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1Cts, startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2Cts, startTime, endTime);
			pmtCDT = getDeadTimeCoinc(pmtCCts,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1Cts.size() + pmt2Cts.size())/(endTime - startTime);
			double rateCoinC = (double)(pmtCCts.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;		
		} else { // "Normal" running
			// Add RDE
			// Deadtime Correction:
			pmt1DT = getDeadTimeSing( pmt1CtsT,startTime, endTime);
			pmt2DT = getDeadTimeSing( pmt2CtsT,startTime, endTime);
			pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
			
			// Pileup Correction
			double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
			double rateCoinC = (double)(coincT.size())/(endTime - startTime);
			Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
			double pmtCPU = puObjC.CalculatePileup(coincT);
			
			pmt1Size = ((double)pmt1CtsT.size()) + pmt1DT + pmt1PU;
			pmt2Size = ((double)pmt2CtsT.size()) + pmt2DT + pmt2PU;
			pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
		}
		
		// Get PMT temp averages (it'll return -1 if no EMS data can be found)
		double temp1 = ems0->getEMSAvg(3311,0,startTime,endTime);
		double temp2 = ems0->getEMSAvg(3311,1,startTime,endTime);
				
		printf("\n BACKGROUNDS --- Run: %d\nheight: %d\nstartTime: %f\nendTime: %f\npmt1cts: %f\npmt2Cts: %f\ncoincCts: %f\nPMT Temps: %f, %f\n",
				runNo, heights[i-1], startTime, endTime, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
		fprintf(outfile,"%d,%d,%f,%f,%f,%f,%f,%f,%f\n",
				runNo, heights[i-1], startTime, endTime, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
		printf("RDE: \npmt1DT: %f\npmt1PU: %f\npmt2DT: %f\npmt2PU: %f\npmtCDT: %f\npmtCPU: %f\n",
				runNo, heights[i-1],pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
		fprintf(outfileRDE,"%d,%d,%f,%f,%f,%f,%f,%f\n",
				runNo, heights[i-1],pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
	}
	fclose(outfile);
	fclose(outfileRDE);
	return;
}

/* Separate out the background for the long run's hold */
void longRunBkg(Run* mcs1, Run* mcs2, EMS* ems0) {
 
	// Figure out when the hold occurs:
	// It begins after cleaning period
	int runNo = mcs1->getRunNo();
	double fillEnd = getFillEnd(mcs1, mcs2);
	double holdStart = mcs2->getTagBitEvt(1<<3, fillEnd-10.0, 1); // Motion of Giant Cleaner
	
	// Check to make sure we have a cleaning period
	if ((holdStart - fillEnd) < 1.0) {
		printf("\nWARNING !!! \nNo cleaning period found, assuming it's 50s !!!\nContinuing anyways...\n");
		holdStart = fillEnd+50.0;
	}
	
	// And it ends when dagger moves.
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	if (dagSteps.front() == dagSteps.back()) {
		printf("Warning !!! \nCannot find dagger dips!!!\n"); 
		return;
	}
	double endTime = dagSteps.front();
	
	// Only want background calibration on long holds
	double startTime = 0.0;
	if ((endTime - holdStart) < 1000.0) {
		printf("Short hold (%f s) !!! Not calculating pre-unload background!!!\n", endTime - holdStart);
		return;
	} else {
		// Counting 70 s before unload starts (This time is 100% clear of background beam)
		startTime = endTime - 70.0;
	}
	
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
	
	
	// Figure out our channels and type of background
	double pmt1Size = 0; // initialize at zero
	double pmt2Size = 0;
	double pmtCSize = 0;
	
	// Rate Dependent components
	double pmt1DT = 0;
	double pmt1PU = 0; // Stays as 0 except for integrated window!
	double pmt2DT = 0;
	double pmt2PU = 0; // Stays as 0 except for integrated window!
	double pmtCDT = 0;
	double pmtCPU = 0;
	
	std::vector<input_t> pmt1CtsT = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1, startTime, endTime](input_t x)->bool{return (x.ch==c1 && x.realtime > startTime && x.realtime < endTime);}
	);
	std::vector<input_t> pmt2CtsT = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2, startTime, endTime](input_t x)->bool{return (x.ch==c2 && x.realtime > startTime && x.realtime < endTime);}
	);
	std::vector<coinc_t> coincT = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[startTime, endTime](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
	);
	// Now we modify our counts based on cMode.
	if ((cMode == 5) || (cMode == 6)) { // Integrated window 
		// Load coincidence info from MCS box
		double coinWindow = dagMCS->getCoincWindow(); // Initial window
		double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
		int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
		
		std::vector<coinc_t> pmt1Cts = getSelfCoincs(pmt1CtsT, coinWindow, PEWindow, PESum);
		std::vector<coinc_t> pmt2Cts = getSelfCoincs(pmt2CtsT, coinWindow, PEWindow, PESum);
		
		// Deadtime Correction:
		pmt1DT = getDeadTimeCoinc(pmt1Cts, startTime, endTime);
		pmt2DT = getDeadTimeCoinc(pmt2Cts, startTime, endTime);
		pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
		
		// Pileup Correction: divided into pmt1 and pmt2 (and combined, duh)
		double ratePmt1  = (double)(pmt1CtsT.size())/(endTime - startTime);
		double rateCoin1 = (double)(pmt1Cts.size()) /(endTime - startTime);
		Pileup puObj1(dagMCS, ratePmt1, rateCoin1, startTime, endTime);
		puObj1.LoadCoincStatsVec(pmt1Cts);
		pmt1PU = puObj1.CalculatePileup(pmt1Cts);
		
		double ratePmt2  = (double)(pmt2CtsT.size())/(endTime - startTime);
		double rateCoin2 = (double)(pmt2Cts.size()) /(endTime - startTime);
		Pileup puObj2(dagMCS, ratePmt2, rateCoin2, startTime, endTime);
		puObj2.LoadCoincStatsVec(pmt2Cts);
		pmt2PU = puObj2.CalculatePileup(pmt2Cts);
		
		double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
		double rateCoinC = (double)(coincT.size())/(endTime - startTime);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
		double pmtCPU = puObjC.CalculatePileup(coincT);
		
		pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
		pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
		pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
	} else if ((cMode == 7) || (cMode == 8)) { // Extra artificial deadtime
		std::vector<input_t> pmt1Cts = imposeDeadtime(pmt1CtsT,500*NANOSECOND);
		std::vector<input_t> pmt2Cts = imposeDeadtime(pmt2CtsT,500*NANOSECOND);
		
		// Add RDE
		// Deadtime Correction:
		pmt1DT = getDeadTimeSing( pmt1Cts, startTime, endTime);
		pmt2DT = getDeadTimeSing( pmt2Cts, startTime, endTime);
		pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
		
		// Pileup Correction
		double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
		double rateCoinC = (double)(coincT.size())/(endTime - startTime);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
		double pmtCPU = puObjC.CalculatePileup(coincT);
		
		pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
		pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
		pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
	} else if ((cMode == 9) || (cMode == 10)) {
		std::vector<input_t> pmt1Cts = removeElectricNoise_sing(pmt1CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<input_t> pmt2Cts = removeElectricNoise_sing(pmt2CtsT,coincT,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<coinc_t> pmtCCts = removeElectricNoise_coinc(coincT,24*NANOSECOND,dagMCS->getPeSum());
		
		// Add RDE
		// Deadtime Correction:
		pmt1DT = getDeadTimeSing( pmt1Cts, startTime, endTime);
		pmt2DT = getDeadTimeSing( pmt2Cts, startTime, endTime);
		pmtCDT = getDeadTimeCoinc(pmtCCts,  startTime, endTime);
		
		// Pileup Correction
		double ratePmtC  = (double)(pmt1Cts.size() + pmt2Cts.size())/(endTime - startTime);
		double rateCoinC = (double)(pmtCCts.size())/(endTime - startTime);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
		double pmtCPU = puObjC.CalculatePileup(coincT);
		
		pmt1Size = ((double)pmt1Cts.size()) + pmt1DT + pmt1PU;
		pmt2Size = ((double)pmt2Cts.size()) + pmt2DT + pmt2PU;
		pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;		
	} else { // "Normal" running
		// Add RDE
		// Deadtime Correction:
		pmt1DT = getDeadTimeSing( pmt1CtsT,startTime, endTime);
		pmt2DT = getDeadTimeSing( pmt2CtsT,startTime, endTime);
		pmtCDT = getDeadTimeCoinc(coincT,  startTime, endTime);
		
		// Pileup Correction
		double ratePmtC  = (double)(pmt1CtsT.size() + pmt2CtsT.size())/(endTime - startTime);
		double rateCoinC = (double)(coincT.size())/(endTime - startTime);
		Pileup puObjC(dagMCS, ratePmtC, rateCoinC, startTime, endTime);
		double pmtCPU = puObjC.CalculatePileup(coincT);
		
		pmt1Size = ((double)pmt1CtsT.size()) + pmt1DT + pmt1PU;
		pmt2Size = ((double)pmt2CtsT.size()) + pmt2DT + pmt2PU;
		pmtCSize = ((double)coincT.size()) + pmtCDT + pmtCPU;
	}
	// Get PMT temp averages (it'll return -1 if no EMS data can be found)
	double temp1 = ems0->getEMSAvg(3311,0,startTime,endTime);
	double temp2 = ems0->getEMSAvg(3311,1,startTime,endTime);
		
	// Output data
	FILE* outfile;
	outfile = fopen("bkgAll.csv", "a");
	//outfile = fopen("longRunBkg.csv", "a");
	printf("Outputting holding backgrounds to longRunBkg.csv !!!\n");
	
	FILE* outfileRDE;
	outfileRDE = fopen("bkgRDE.csv","a");
	printf("Outputting RDE data to bkgRDE.csv!\n");
	printf("\n BACKGROUNDS --- Run: %d\nheight: %d\nstartTime: %f\nendTime: %f\npmt1cts: %f\npmt2Cts: %f\ncoincCts: %f\nPMT Temps: %f, %f\n",
		runNo, 490, startTime, endTime, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
	fprintf(outfile,"%d,%d,%f,%f,%f,%f,%f,%f,%f\n",
		runNo, 490, startTime, endTime, pmt1Size, pmt2Size, pmtCSize, temp1, temp2);
	printf("RDE: \npmt1DT: %f\npmt1PU: %f\npmt2DT: %f\npmt2PU: %f\npmtCDT: %f\npmtCPU: %f\n",
			runNo, 490,pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
	fprintf(outfileRDE,"%d,%d,%f,%f,%f,%f,%f,%f\n",
			runNo, 490,pmt1DT,pmt1PU,pmt2DT,pmt2PU,pmtCDT,pmtCPU);
	fclose(outfile);
	fclose(outfileRDE);
	return;
}

/* This is the "background" during normal production running.
 * "Background" periods are 20s during hold, 50s at end of dagger, and
 * 40s after the last cat door event.  */
void productionBkg(Run* mcs1, Run* mcs2, EMS* ems0) {
 
	// Figure out when the hold occurs:
	// It begins after cleaning period
	// Load tagbit data:
	int runNo = mcs1->getRunNo();
	double fillEnd = getFillEnd(mcs1, mcs2);
			
	// Make sure the run actually loaded...
	// This is just to keep the same runs in "production" as here.
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2, fillEnd);
	if(beamHits.empty() || (beamHits.end() <= beamHits.begin())) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", runNo);
		return;
	}
	// Find the dagger steps
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
	if(dagSteps.empty() || (dagSteps.end() <= dagSteps.begin())) { 
		fprintf(stderr, "Skipping run %d for no Dagger Movement\n", runNo);
		return;
	}
	
	// Now let's look at timings.
	double holdEnd = dagSteps.front(); // Hold ends when dagger moves.
	double holdStart = holdEnd - 20.0; // Assume 20s (the shortest) hold
	double bkgEnd   = *(dagSteps.end()-1); // End of the background period
	double bkgStart = *(dagSteps.end()-1) - 40; // Start of background
	double cntEnd   = bkgStart  - 10; // End of unload	
	if (dagSteps.size() > 2) { // Avoid segfaults, get slightly more accurate guesses
		bkgStart = mcs1->getTagBitEvt(1<<4, *(dagSteps.end()-2), 1) + 5.0;// Cat Door movement + 5
		cntEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1); // Trap door movement
	}
	double cntStart = cntEnd - 50; // Last 50 seconds of unload

	// Double check these are all in order:
	if ((holdStart <= fillEnd) || (cntStart <= holdEnd) || (cntEnd <= cntStart) ||
		(bkgStart <= cntEnd) || (bkgEnd <= bkgStart)) {
		fprintf(stderr, "Skipping run %d for background timing glitch\n", runNo);
		return;
	}
	
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
		
	// Singles counts (for variable integration windows)
	// Coincidences are always type coinc_t so there's not an issue there.
	double hold1Size = 0; // initialize at zero
	double hold2Size = 0;
	double holdCSize = 0;
	double cnt1Size  = 0;
	double cnt2Size  = 0;
	double cntCSize  = 0;
	double bkg1Size  = 0;
	double bkg2Size  = 0;
	double bkgCSize  = 0;
	
	std::vector<input_t> hold1Tmp = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1, holdStart, holdEnd](input_t x)->bool{return (x.ch==c1 && x.realtime > holdStart && x.realtime < holdEnd);}
	);
	std::vector<input_t> hold2Tmp = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2, holdStart, holdEnd](input_t x)->bool{return (x.ch==c2 && x.realtime > holdStart && x.realtime < holdEnd);}
	);
	std::vector<coinc_t> holdCTmp = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[holdStart, holdEnd](coinc_t x)->bool{return x.realtime > holdStart && x.realtime < holdEnd;}
	);
	std::vector<input_t> cnt1Tmp = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1, cntStart, cntEnd](input_t x)->bool{return (x.ch==c1 && x.realtime > cntStart && x.realtime < cntEnd);}
	);
	std::vector<input_t> cnt2Tmp = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2, cntStart, cntEnd](input_t x)->bool{return (x.ch==c2 && x.realtime > cntStart && x.realtime < cntEnd);}
	);
	std::vector<coinc_t> cntCTmp = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[cntStart, cntEnd](coinc_t x)->bool{return x.realtime > cntStart && x.realtime < cntEnd;}
	);
	std::vector<input_t> bkg1Tmp = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1, bkgStart, bkgEnd](input_t x)->bool{return (x.ch==c1 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	std::vector<input_t> bkg2Tmp = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2, bkgStart, bkgEnd](input_t x)->bool{return (x.ch==c2 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	std::vector<coinc_t> bkgCTmp = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[bkgStart, bkgEnd](coinc_t x)->bool{return x.realtime > bkgStart && x.realtime < bkgEnd;}
	);
			
	// Now let's load all these singles vectors.
	if (cMode < 5) { // Normal Singles running
		
		hold1Size = hold1Tmp.size(); // Due to compiler issues, put these inside the if
		hold2Size = hold2Tmp.size();
		holdCSize = holdCTmp.size();
		cnt1Size = cnt1Tmp.size();
		cnt2Size = cnt2Tmp.size();
		cntCSize = cntCTmp.size();
		bkg1Size = bkg1Tmp.size();
		bkg2Size = bkg2Tmp.size();
		bkgCSize = bkgCTmp.size();
	} else if ((cMode == 5) || (cMode == 6)) { // integrated window 
		// Load coincidence info for PMT
		double coinWindow = dagMCS->getCoincWindow(); // Initial window
		double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
		int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
		// And convert to coincidence
		std::vector<coinc_t> hold1Cts = getSelfCoincs(hold1Tmp, coinWindow, PEWindow, PESum);
		std::vector<coinc_t> hold2Cts = getSelfCoincs(hold2Tmp, coinWindow, PEWindow, PESum);
		std::vector<coinc_t> cnt1Cts = getSelfCoincs(cnt1Tmp, coinWindow, PEWindow, PESum);
		std::vector<coinc_t> cnt2Cts = getSelfCoincs(cnt2Tmp, coinWindow, PEWindow, PESum);
		std::vector<coinc_t> bkg1Cts = getSelfCoincs(bkg1Tmp, coinWindow, PEWindow, PESum);
		std::vector<coinc_t> bkg2Cts = getSelfCoincs(bkg2Tmp, coinWindow, PEWindow, PESum);
		
		hold1Size = hold1Cts.size(); // Due to compiler issues, put these inside the if
		hold2Size = hold2Cts.size();
		holdCSize = holdCTmp.size();
		cnt1Size = cnt1Cts.size();
		cnt2Size = cnt2Cts.size();
		cntCSize = cntCTmp.size();
		bkg1Size = bkg1Cts.size();
		bkg2Size = bkg2Cts.size();
		bkgCSize = bkgCTmp.size();
	} else if ((cMode == 7) || (cMode == 8)) { // Extra artificial deadtime
		std::vector<input_t> hold1Cts = imposeDeadtime(hold1Tmp,500*NANOSECOND);
		std::vector<input_t> hold2Cts = imposeDeadtime(hold2Tmp,500*NANOSECOND);
		std::vector<input_t> cnt1Cts = imposeDeadtime(cnt1Tmp,500*NANOSECOND);
		std::vector<input_t> cnt2Cts = imposeDeadtime(cnt2Tmp,500*NANOSECOND);
		std::vector<input_t> bkg1Cts = imposeDeadtime(bkg1Tmp,500*NANOSECOND);
		std::vector<input_t> bkg2Cts = imposeDeadtime(bkg2Tmp,500*NANOSECOND);

		hold1Size = hold1Cts.size(); // Due to compiler issues, put these inside the if
		hold2Size = hold2Cts.size();
		holdCSize = holdCTmp.size();
		cnt1Size = cnt1Cts.size();
		cnt2Size = cnt2Cts.size();
		cntCSize = cntCTmp.size();
		bkg1Size = bkg1Cts.size();
		bkg2Size = bkg2Cts.size();
		bkgCSize = bkgCTmp.size();

	} else if ((cMode == 9) || (cMode == 10)) {
		std::vector<input_t> hold1Cts = removeElectricNoise_sing(hold1Tmp,holdCTmp,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<input_t> hold2Cts = removeElectricNoise_sing(hold2Tmp,holdCTmp,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<coinc_t> holdCCts = removeElectricNoise_coinc(holdCTmp,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<input_t> cnt1Cts = removeElectricNoise_sing(cnt1Tmp,cntCTmp,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<input_t> cnt2Cts = removeElectricNoise_sing(cnt2Tmp,cntCTmp,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<coinc_t> cntCCts = removeElectricNoise_coinc(cntCTmp,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<input_t> bkg1Cts = removeElectricNoise_sing(bkg1Tmp,bkgCTmp,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<input_t> bkg2Cts = removeElectricNoise_sing(bkg2Tmp,bkgCTmp,24*NANOSECOND,dagMCS->getPeSum());
		std::vector<coinc_t> bkgCCts = removeElectricNoise_coinc(bkgCTmp,24*NANOSECOND,dagMCS->getPeSum());
		
		hold1Size = hold1Cts.size(); // Due to compiler issues, put these inside the if
		hold2Size = hold2Cts.size();
		holdCSize = holdCCts.size();
		cnt1Size = cnt1Cts.size();
		cnt2Size = cnt2Cts.size();
		cntCSize = cntCCts.size();
		bkg1Size = bkg1Cts.size();
		bkg2Size = bkg2Cts.size();
		bkgCSize = bkgCCts.size();
	}
		
	// Output 2 .csv files that we will use in NormAnalyzer-Coinc.py
	// Start by outputting monitor and backgrounds
	FILE* outfileBkg;
	outfileBkg = fopen("production-bkg.csv", "a");
	printf("Outputting background data to production-bkg.csv!\n");
	printf("Monitor Data -----Run: %d\ntdTime: %f\nholdStart: %f\nholdEnd: %f\nhold1: %f\nhold2: %f\nholdC: %f\ncntStart: %f\ncntEnd: %f\ncnt1: %f\ncnt2: %f\ncntC: %f\nbkgStart: %f\nbkgEnd: %f\nbkg1: %f\nbkg2: %f\nbkgC: %f\n",
			runNo, fillEnd, 
			holdStart,holdEnd,hold1Size,hold2Size,holdCSize,
			cntStart,cntEnd,cnt1Size,cnt2Size,cntCSize,
			bkgStart,bkgEnd,bkg1Size,bkg2Size,bkgCSize);
	fprintf(outfileBkg,"%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
			runNo, fillEnd, 
			holdStart,holdEnd,hold1Size,hold2Size,holdCSize,
			cntStart,cntEnd,cnt1Size,cnt2Size,cntCSize,
			bkgStart,bkgEnd,bkg1Size,bkg2Size,bkgCSize);
	fclose(outfileBkg);
}
