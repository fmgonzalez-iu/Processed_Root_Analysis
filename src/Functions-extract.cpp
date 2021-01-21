#include "../inc/Functions.hpp"

/* define constants we need for later */
#define NANOSECOND .000000001
#define WINDOW 40000
#define bkgMov50ns8pe 0.10666
#define bkgMov50ns1000ns6pe 0.3
#define synthbkg_50_500_2 0.0
#define TAUN 877.7

/*----------------------------------------------------------------------
 *		Author: F. Gonzalez
 * 
 * 		These are all the functions that extract values from ROOT's 
 * 		processed_output_xxxxx.root files into the terminal/.csv for
 * 		analysis. 
 * 
 *--------------------------------------------------------------------*/

/* The tagbitTagger (TM) searches through data and determines what runs
 * have "good" tagbits. */
int tagbitTagger(Run* mcs1, Run* mcs2, bool verbose) {
	/*Tagbit table (from 2017-2018): 
	MCS_0:
	* 0) Small Cleaner
	* 1) N/A
	* 2) Gate Valve
	* 3) Trap Door
	* 4) Cat Door
	MCS_1:
	* 0) Dagger Movement
	* 1) Beam On
	* 2) Giant Cleaner
	* 3) N/A
	* 4) N/A*/
	
	// Load run number
	int runNo = mcs1->getRunNo();
	
	// 2 second timing slop (arbitrary)
	double slop = 2.0;
	
	// These events should happen once per run
	double acUp    = mcs1->getTagBitEvt(1<<0, slop, 1);
	double gvClose = mcs1->getTagBitEvt(1<<2, slop, 1);
	double tdClose = mcs1->getTagBitEvt(1<<3, slop, 1);
	double cdClose = mcs1->getTagBitEvt(1<<4, slop, 1);
	double gcUp    = mcs2->getTagBitEvt(1<<2, slop, 1); // GC doesn't really tag
	
	// These events might not happen every run, but could.
	double acDown  = mcs1->getTagBitEvt(1<<0, acUp + slop, 0);
	double tdOpen  = mcs1->getTagBitEvt(1<<3, tdClose + slop, 1);
	double cdOpen  = mcs1->getTagBitEvt(1<<4, cdClose + slop, 1);
	double gcDown  = mcs2->getTagBitEvt(1<<2, gcUp + slop, 1); // GC doesn't really tag
	
	// H-GX and dagger have multiple (variable) status changes.
	double loopStart = 0.0;
	double loopStop  = 0.0;
	
	// Find all times the dagger moves during the run
	std::vector<double> dagMoves;
	std::vector<double> dagStops;
	loopStart = mcs2->getTagBitEvt(1<<0, slop, 1);
	do {
		dagMoves.push_back(loopStart);
		loopStop = mcs2->getTagBitEvt(1<<0, dagMoves.back() + slop, 0);
		dagStops.push_back(loopStop);
		loopStart= mcs2->getTagBitEvt(1<<0, loopStop + slop, 1);
	} while(loopStart != -1.0 && loopStop != -1.0);
            
    // Look for all HGX pulses in the DAQ
    std::vector<double> HGXPulses;
    loopStart = 0.0;
    loopStop  = 0.0;
	while(loopStart < dagMoves.back() && loopStart >= 0.0) {
		HGXPulses.push_back(loopStart);
		loopStart = mcs2->getTagBitEvt(1<<1, loopStart + slop, 1);
	}

	if (verbose) {
	
		// Verbose will just print out the times (for debugging purposes)
		printf("\nOUTPUTTING TIMINGS FOR RUN %d:\n\n", runNo);
		printf("   Active Cleaner raises at %f and lowers again at %f\n", acUp, acDown);
		printf("   Gate Valve closes at %f\n", gvClose);
		printf("   Trap Door closes at %f and opens at %f\n", tdClose, tdOpen);
		printf("   Cat Door closes at %f and opens at %f\n", cdClose, cdOpen);
		printf("   Giant Cleaner raises at %f and lowers again at %f\n", gcUp, gcDown);
		printf("   Dagger initially raises at %f\n", dagMoves[0]);
		printf("   This run has %lu unload steps:\n", dagMoves.size());
		if(dagMoves.size() > 0) {
			for (auto dIt=dagMoves.begin(); dIt<dagMoves.end();dIt++){
				printf("      Dagger movement at %f\n",(*dIt));
			}
		}
		if(dagStops.size() > 0) {
			for (auto dIt=dagStops.begin();dIt<dagStops.end();dIt++){
				printf("      Dagger stops at %f\n",(*dIt));
			}
		}
	}
	
	// Find average HGX spacing:
	int fillPulses = 0;
	double fillSpacing = 0.0;
	if (HGXPulses.size() > 0){
		for (auto hgx = HGXPulses.begin()+1; hgx < HGXPulses.end(); hgx++) {
			// Calculating spacing only while GV is open
			if ((*hgx < gvClose) || (gvClose < 0.0)) { 
				fillSpacing += (*hgx) - *(hgx-1);
				fillPulses += 1;
			}		
		}
	}
	// normalize pulse spacing
	if (fillPulses > 0) {
		fillSpacing = fillSpacing/((double)fillPulses);
	} else {
		fillSpacing = -1;
	}
	
	if (verbose && HGXPulses.size() > 0) {
		printf("   There are %d H-GX hits, starting at %f and ending at %f.\n", fillPulses, HGXPulses.front(),HGXPulses[fillPulses]);	
		printf("    HGX spacing is %f.\n", fillSpacing); 
	}
	
	if (HGXPulses.size() < 10) {
		printf("Not enough H-GX pulses! Probably no beam on for this run!\n");
		return -1;
	}
	
	// number of dagger unload dips
	int dips = dagMoves.size() - 2;
	if (dips != 3) {
		printf("Not normal production unload! Might be 1-dip or 9-dip (or a bug...)!\n");
	}

	double cleanT = acUp - gvClose;
	if (cleanT <= 0) {
		cleanT = 50.0; // Error, so defaulting to 50
	}
	double holdT  = dagMoves[1]-(gvClose+cleanT);
	if (holdT <= slop) {
		holdT = -1.0; // Error
	}
	double bkgT   = dagMoves.back() - tdOpen;
	if (bkgT <= slop) {
		bkgT = -1.0; // Error
	}
	FILE* outfile;
	outfile = fopen("tagTiming.csv", "a");
	printf("Outputting data to tagTiming.csv!\n");

	printf("Run: %d \nFill: %f \nClean: %f \nHold: %f \nnDips: %d \nH-GX Spacing: %f \nH-GX Pulses: %d \nBackground: %f\n",runNo,gvClose,cleanT,holdT,dips,fillSpacing,fillPulses,bkgT);
	fprintf(outfile,"%d,%f,%f,%f,%d,%f,%d,%f",runNo,gvClose,cleanT,holdT,dips,fillSpacing,fillPulses,bkgT);
	fclose(outfile);
	return 0;
}

void holdingPressure(Run* mcs1, Run* mcs2, EMS* ems0) {
	
	// Load tagbit data: FillEnd, H-GX, and Dagger Movement
	int runNo = mcs1->getRunNo();
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
	
	double countEnd;
	if (dagSteps.size() > 2) { 
		countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
	} else {
		countEnd = *(dagSteps.end()-1);
	}
	

	double cleaner  = fillEnd;
	double counting = countEnd;
	FILE* outfile = fopen("mean_trap_pressures_wClean.csv","a");
	// Get PMT temp averages (it'll return -1 if no EMS data can be found)
	double presN = ems0->getEMSAvg(1200,50,cleaner,counting);
	double minN  = ems0->getEMSMin(1200,50,cleaner,counting);
	double maxN  = ems0->getEMSMax(1200,50,cleaner,counting);
	double presS = ems0->getEMSAvg(1200,52,cleaner,counting);
	double minS  = ems0->getEMSMin(1200,52,cleaner,counting);
	double maxS  = ems0->getEMSMax(1200,52,cleaner,counting);
	
	printf("Data for run %d:\nt1: %f\nt2: %f\npressureN: %e (%e to %e) \npressureS: %e (%e to %e)\n",runNo,cleaner,counting,presN,minN,maxN,presS,minS,maxS);
	fprintf(outfile,"%d,%f,%f,%e,%e,%e,%e,%e,%e\n",runNo,cleaner,counting,presN,minN,maxN,presS,minS,maxS);
	
	fclose(outfile);
	
}


int coincidenceCompare(Run* mcs1, Run* mcs2) {
	// For running comparisons between 
	
	Run* dagMCS;
	int cMode = mcs1->getCoincMode();
	// Determine which dagger pair we want to use
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	double start = -1.0;
	double end = 999999.0;
	
	std::vector<coinc_t> coincs = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[start, end](coinc_t x)->bool{return (x.realtime > start && x.realtime < end);}
	);
		
	FILE* outfile;
	//outfile = fopen("testList_FMG.csv","a");
	outfile = fopen("coincList_FMG.csv_long","w");
	//unsigned long begin = 30824423355; // 
	//unsigned long begin = 79202002712; //
	unsigned long begin = 357430180467; //
	//unsigned long length = 2583; // 
	//unsigned long length = 2031; //
	unsigned long length = 1958; //
	std::vector<input_t> input = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1,c2,begin, length](input_t x)->bool{return ((x.ch == c1 || x.ch == c2) && x.time >= begin && x.time <= begin+length);}
	);
	/*for (auto cIt=input.begin();cIt<input.end();cIt++){
		fprintf(outfile,"%llu,%d\n",(*cIt).time,(*cIt).ch);
	}*/
	for (auto cIt=coincs.begin(); cIt<coincs.end();cIt++) {
		unsigned long len;
		len = (unsigned long)((*cIt).length / (0.8*NANOSECOND) + 0.5);
		fprintf(outfile,"%llu,%d,%d,%llu\n",(*cIt).time,(*cIt).pmt1,(*cIt).pmt2,len);
	}
	return (int)coincs.size();
}

//----------------------------------------------------------------------
// Important Code!!!
//----------------------------------------------------------------------
/* normNByDip simultaneously does singles/coincidence runs.
 * Basically, it's a pseudo-optimized 
 * normNByDipDet->normNByDipCoinc/normNByDipSing. */
void normNByDip(Run* mcs1, Run* mcs2, EMS* ems0, const char* filename) {
	
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

	// load time constant data
	std::vector<int> rawDet;
	std::vector<double> kappa; 
	std::vector<double> kappaE;
	FILE* tc_file = fopen(filename, "r");
	if (tc_file != NULL) {
		char line[48];
		while (fgets(line,48,tc_file)) {
			// raw detectors coming from python script as 1-10. 
			const char* minRun = strtok(line,",");
			const char* maxRun = strtok(NULL,",");
			if ((std::stoi(minRun) <= runNo) && (std::stoi(maxRun) > runNo)) {
				const char* detF = strtok(NULL,",");
				rawDet.push_back(std::stoi(detF));
				kappa.push_back(std::stod(strtok(NULL,",")));
				kappaE.push_back(std::stod(strtok(NULL,",")));
			}
		}
		fclose(tc_file);
	}
	if (rawDet.size() == 0) { // If there's a problem in loading tc_file assume {N, 0.0,0.0}
		//printf("No TC file found, assuming 70 s fit!\n");
		for (int i = 1; i < 11; i++) {
			rawDet.push_back(i);
			kappa.push_back(0.0); // No tc_file means no weighting
			kappaE.push_back(0.0);
		}
	}
	std::vector<double> mBkg;
	std::vector<double> mBkgE;
	const char* monBkgLoc  = std::getenv("MONBKG_LOC"); 
	FILE* mb_file = fopen(monBkgLoc,"r");
	if (mb_file != NULL) {
		char line[48];
		while (fgets(line,48,mb_file)) {
			// Monitor Backgrounds coming from python script as 1-10.
			const char* minRun = strtok(line,",");
			const char* maxRun = strtok(NULL,",");
			if ((std::stoi(minRun) <= runNo) && (std::stoi(maxRun) > runNo)) {
				const char* detF = strtok(NULL,",");
				mBkg.push_back(std::stod(strtok(NULL,",")));
				mBkgE.push_back(std::stod(strtok(NULL,",")));
			}
		}
		fclose(mb_file);
	}	
	if (rawDet.size() != mBkg.size()) {
		mBkg.clear();
		mBkgE.clear();
		for (size_t i = 0; i < rawDet.size(); i++){
			mBkg.push_back(0.0);
			mBkgE.push_back(0.0);
		}
	}
		
	// Load a vector of weighted monitors.
	Run* dagMCS;
	std::vector<measurement> weightMon;
	std::vector<input_t> detCts;
	for (int ii = 0; ii < rawDet.size(); ii++) {
		int detCh;
		// need to convert detectors 1-10 to the right MCS box
		if (ii < 5) {
			dagMCS = mcs1;
			detCh = rawDet[ii] + dagMCS->getMCSOff();
		} else {
			dagMCS = mcs2;
			detCh = (rawDet[ii]-5) + dagMCS->getMCSOff();
		}
		detCts = dagMCS->getCounts(
			//[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
			[](input_t x)->input_t{return x;},
			[detCh, fillEnd](input_t x)->bool{return x.ch == detCh && x.realtime < fillEnd;}
		);
	
		//weightMon.push_back(expWeightMonVectDT(detCts, kappa[ii],16*NANOSECOND,fillEnd));
		weightMon.push_back(expWeightMonVectBkg(detCts, kappa[ii],16*NANOSECOND,mBkg[ii],fillEnd));
	}
	printf("Loaded Monitor Counts\n");
	// Find background times
	// During production running, the last 50s should have TD open, dagger down bkgs
	double lastDip;
	double countEnd;
	double bkgStart;
	double bkgEnd = *(dagSteps.end()-1);
	if (dagSteps.size() > 2) { 
		lastDip = *(dagSteps.end() - 2);
		countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
		bkgStart = mcs1->getTagBitEvt(1<<4, *(dagSteps.end()-2), 1) + 5.0;// 1<<4 is Cat Door, add extra 5 seconds to be safe
	} else {
		bkgStart = *(dagSteps.end()-1) - 40;
		fprintf(stderr, "\nUsing LAST 40s for background!!!!!!!\n");
	}
	
	if (bkgStart < 5.0 || bkgStart >= bkgEnd) { // tag = -1.0 for errors (we've added 5s for stability.)
		fprintf(stderr, "Skipping run %d for background timing glitch\n", runNo);
		return;
	}
	
	// Determine which dagger pair we want to use
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	//bkgStart = *dagSteps.begin();
	// Load background counts -- first temp
	std::vector<input_t> bkgDag1t = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c1 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	std::vector<input_t> bkgDag2t = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c2 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	printf("Loaded singles background\n");
	std::vector<coinc_t> bkgDagCt = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[bkgStart, bkgEnd](coinc_t x)->bool{return (x.realtime > bkgStart && x.realtime < bkgEnd);}
	);			
	printf("Loaded coinces background\n");
	// Now do operations on our temp counts
	//std::vector<input_t> bkgDag1;
	//std::vector<input_t> bkgDag2;
	//std::vector<coinc_t> bkgDagC;
	
	double bkgDag1Size;
	double bkgDag2Size;
	double bkgDagCSize;
	if ((cMode == 5) || (cMode == 6)) { // Integrated window
		double coinWindow = dagMCS->getCoincWindow(); // Initial window
		double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
		int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
		// Need to RDE these...
		std::vector<coinc_t> bkgDag1 = getSelfCoincs(bkgDag1t, coinWindow, PEWindow, PESum);
		bkgDag1Size = (double)bkgDag1.size();
		std::vector<coinc_t> bkgDag2 = getSelfCoincs(bkgDag2t, coinWindow, PEWindow, PESum);
		bkgDag2Size = (double)bkgDag2.size();
		std::vector<coinc_t> bkgDagC = bkgDagCt;
		bkgDagCSize = (double)bkgDagC.size();
	} else if ((cMode == 7) || (cMode == 8)) {
		std::vector<input_t> bkgDag1 = imposeDeadtime(bkgDag1t,500*NANOSECOND);
		bkgDag1Size = (double)bkgDag1.size();
		std::vector<input_t> bkgDag2 = imposeDeadtime(bkgDag2t,500*NANOSECOND);
		bkgDag2Size = (double)bkgDag2.size();
		std::vector<coinc_t> bkgDagC = bkgDagCt;
		bkgDagCSize = (double)bkgDagC.size();
	} else if ((cMode == 9) || (cMode == 10)) {
		std::vector<input_t> bkgDag1 = removeElectricNoise_sing(bkgDag1t,bkgDagCt,24*NANOSECOND,dagMCS->getPeSum());
		bkgDag1Size = (double)bkgDag1.size();
		std::vector<input_t> bkgDag2 = removeElectricNoise_sing(bkgDag2t,bkgDagCt,24*NANOSECOND,dagMCS->getPeSum());
		bkgDag2Size = (double)bkgDag2.size();
		std::vector<coinc_t> bkgDagC = removeElectricNoise_coinc(bkgDagCt,24*NANOSECOND,dagMCS->getPeSum());
		bkgDagCSize = (double)bkgDagC.size();
		
		printf("Removed Noise\n");
		//bkgDag1Size += (unsigned long)ceil(getDeadTimeSing(bkgDag1, bkgStart, bkgEnd));
		//bkgDag2Size += (unsigned long)ceil(getDeadTimeSing(bkgDag2, bkgStart, bkgEnd));
		//bkgDagCSize += (unsigned long)ceil(getDeadTimeCoinc(bkgDagC, bkgStart, bkgEnd));	
		bkgDag1Size += getDeadTimeSing(bkgDag1, bkgStart, bkgEnd);
		bkgDag2Size += getDeadTimeSing(bkgDag2, bkgStart, bkgEnd);
		bkgDagCSize += getDeadTimeCoinc(bkgDagC, bkgStart, bkgEnd);	
	} else {
		bkgDag1Size = (double)bkgDag1t.size();
		bkgDag2Size = (double)bkgDag2t.size();
		bkgDagCSize = (double)bkgDagCt.size();
	}
	//------------------------------------------------------------------
	// Create a pileup object for accounting for pileup
	// Go across dip 3 for now
	double bkgRaw = (double)(bkgDag1Size+bkgDag2Size)/(bkgEnd-bkgStart);
	double bkgCoinc = (double)(bkgDagCSize)/(bkgEnd-bkgStart);
	
	printf("Loaded Backgrounds, generating PU\n");
	Pileup puObj(dagMCS, bkgRaw,bkgCoinc,lastDip, countEnd);
	// Should probably make this flexible, but ehhh
	// If we're filtering, need to load singles/coinc data here.
	double tau = puObj.GetTau();
	printf("Tau = %e\n", tau);
	double pileup = puObj.CalculatePileup(bkgDagCt);
	printf("PU = %f\n",pileup);
	printf("amp = %f, dt = %f\n",puObj.GetAmplitude(),tau);
	//------------------------------------------------------------------
	// Output 2 .csv files that we will use in NormAnalyzer-Coinc.py
	// Start by outputting monitor and backgrounds
	FILE* outfileDet;
	outfileDet = fopen("normNByDip-det.csv", "a");
	printf("Outputting monitor/background data to normNByDip-det.csv!\n");
	printf("\nMonitor Data -----Run: %d\ntdTime: %f\nholdEnd: %f\nbkgDag1: %f\nbkgDag2: %f\nbkgDagC: %f\nbkgStart: %f\nbkgEnd: %f\n",
			runNo, fillEnd, *dagSteps.begin(), bkgDag1Size,bkgDag2Size,bkgDagCSize+pileup, bkgStart, bkgEnd);
	fprintf(outfileDet,"%d,%f,%f,%f,%f,%f,%f,%f",
			runNo, fillEnd, *dagSteps.begin(), bkgDag1Size,bkgDag2Size,bkgDagCSize+pileup, bkgStart, bkgEnd);			
	for (int jt = 0; jt < weightMon.size(); jt++ ) {
		printf("Monitor %d: %f (%f)\n",jt,weightMon[jt].val,weightMon[jt].err);
		fprintf(outfileDet,",%f,%f",weightMon[jt].val,weightMon[jt].err);
	}
	fprintf(outfileDet,"\n");
	fclose(outfileDet);
	
	// Now output dagger detector data
	FILE* outfileSing;
	outfileSing = fopen("normNByDipSing-dag.csv","a");
	FILE* outfileCoin;
	outfileCoin = fopen("normNByDipCoinc-dag.csv","a");
	
	// Separating dips into 10s slices (arb.)
	int numSteps = 0;
	double sliceTime = 10.0;
	size_t sC = 0; // Counter for dagger stops
	
	
	// TEMPORARY--------------------------------------------------------
	FILE* outfilePU;
	outfilePU = fopen("pileup.csv","a");
	
	double pileupTot = 0.;
	// Loop through each dagger step and find when photon events hit
	for(auto stepIt = dagSteps.begin()+1; stepIt < dagSteps.end(); stepIt++) {
		
		// Time for this step
		double startTime = *(stepIt-1);
		double endTime = *(stepIt);
				
		// Cut off the background period in the last 50s.
		if(stepIt == dagSteps.end()-1) { 
			endTime = countEnd;
		}
		
		if (startTime >= endTime) {
			printf("Stopping due to timing bug!\n");
			return;
		}
		
		// Divide each step into 10 second slices
		double tempSta = startTime;
		double tempEnd = startTime + sliceTime;
		//for (int jj = 0; jj < (int)((endTime - startTime) / sliceTime) + 1; jj++) {
		for (int jj = 0; jj < (int)((endTime - startTime) / sliceTime) + 2; jj++) {
			
			// Want to add a temporary step for the first (moving) time.
			if (jj == 0) {
				tempEnd = dagStop.at(sC);
			}
			// Make sure we don't go too far
			if (tempEnd > (endTime - 0.2)) {
				tempEnd = endTime;
			}
			// 0.2s timing slop
			if ((tempSta + 0.2) > tempEnd) { 
				break;
			}
			
			// Load counts (singles + coincidence)
			std::vector<input_t> cts1t = dagMCS->getCounts(
				[](input_t x)->input_t{return x;},
				[c1,tempSta,tempEnd](input_t x)->bool{return (x.ch == c1 && x.realtime > tempSta && x.realtime < tempEnd);}
			);
			std::vector<input_t> cts2t = dagMCS->getCounts(
				[](input_t x)->input_t{return x;},
				[c2,tempSta,tempEnd](input_t x)->bool{return (x.ch == c2 && x.realtime > tempSta && x.realtime < tempEnd);}
			);
			std::vector<coinc_t> ctsCt = dagMCS->getCoincCounts(
				[](coinc_t x)->coinc_t{return x;},
				[tempSta, tempEnd](coinc_t x)->bool{return (x.realtime > tempSta && x.realtime < tempEnd);}
			);
			
			// And do whatever filtering we need to do
			//std::vector<input_t> cts1;
			//std::vector<input_t> cts2;
			std::vector<coinc_t> ctsC; // Always the same format
			
			// Load variables for output -- internals differ for cModes
			unsigned long cts1Size; // Size of counts vector
			unsigned long cts2Size;
			unsigned long ctsCSize;
			
			double corr1; // Size of RDE missed counts
			double corr2;
			//double corrC;
			
			double stepMean1; // Mean arrival time steps
			double stepMean2;
			//double stepMeanC;
			
			double coin1 = 0.; // Coincidence/Singles conversion
			double coin2 = 0.;
			
			if ((cMode == 5) || (cMode == 6)) { // Integrated window
				double coinWindow = dagMCS->getCoincWindow(); // Initial window
				double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
				int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
		
				std::vector<coinc_t> cts1 = getSelfCoincs(cts1t, coinWindow, PEWindow, PESum);
				std::vector<coinc_t> cts2 = getSelfCoincs(cts2t, coinWindow, PEWindow, PESum);
				ctsC = ctsCt; // For now ctsC is previously declared...
				
				cts1Size = cts1.size();
				cts2Size = cts2.size();
				ctsCSize = ctsC.size();
				// Need to do different 
				stepMean1 = cts1.size() > 0 ? std::accumulate(cts1.begin(), cts1.end(), 0.0, [](double m, coinc_t x)->double{return m + x.realtime;}) / (double)cts1.size() : 0.0;
				stepMean2 = cts2.size() > 0 ? std::accumulate(cts2.begin(), cts2.end(), 0.0, [](double m, coinc_t x)->double{return m + x.realtime;}) / (double)cts2.size() : 0.0;
			
				corr1 = getDeadTimeCoinc( cts1, tempSta, tempEnd);
				corr2 = getDeadTimeCoinc( cts2, tempSta, tempEnd);
				
				if (ctsC.size() > 0) { // divide by zero
					coin1 = ((double)cts1Size)/((double)ctsCSize);
					coin2 = ((double)cts2Size)/((double)ctsCSize);
				}
			} else { // cts1 as format input_t
				std::vector<input_t> cts1;
				std::vector<input_t> cts2;
				if ((cMode == 7) || (cMode == 8)) {
					cts1 = imposeDeadtime(cts1t,500*NANOSECOND);
					cts2 = imposeDeadtime(cts2t,500*NANOSECOND);
					ctsC = ctsCt;
					
				} else if ((cMode == 9) || (cMode == 10)) {
					//printf("removing electric noise, time %f,%f\n",tempSta,tempEnd);
					cts1 = removeElectricNoise_sing(cts1t,ctsCt,24*NANOSECOND,dagMCS->getPeSum());
					cts2 = removeElectricNoise_sing(cts2t,ctsCt,24*NANOSECOND,dagMCS->getPeSum());
					ctsC = removeElectricNoise_coinc(ctsCt,24*NANOSECOND,dagMCS->getPeSum());
					
				} else {
					cts1 = cts1t;
					cts2 = cts2t;
					ctsC = ctsCt;
				}
				cts1Size = cts1.size();
				cts2Size = cts2.size();
				ctsCSize = ctsC.size();
				
				stepMean1 = cts1.size() > 0 ? std::accumulate(cts1.begin(), cts1.end(), 0.0, [](double m, input_t x)->double{return m + x.realtime;}) / (double)cts1.size() : 0.0;
				stepMean2 = cts2.size() > 0 ? std::accumulate(cts2.begin(), cts2.end(), 0.0, [](double m, input_t x)->double{return m + x.realtime;}) / (double)cts2.size() : 0.0;
				corr1 = getDeadTimeSing( cts1, tempSta, tempEnd);
				corr2 = getDeadTimeSing( cts2, tempSta, tempEnd);
				
				if (ctsC.size() > 0) { // divide by zero
					coin1 = 0.0;
					coin2 = 0.0;
					for (auto cIt = ctsC.begin(); cIt < ctsC.end(); cIt++) {
						coin1 += (double)cIt->pmt1;
						coin2 += (double)cIt->pmt2;
					}
					coin1 = coin1 / (double)ctsC.size();
					coin2 = coin2 / (double)ctsC.size();
				}
			}			
			// Find Mean Arrival Time
			//double stepMean1 = cts1.size() > 0 ? std::accumulate(cts1.begin(), cts1.end(), 0.0, [](double m, input_t x)->double{return m + x.realtime;}) / (double)cts1.size() : 0.0;
			//double stepMean2 = cts2.size() > 0 ? std::accumulate(cts2.begin(), cts2.end(), 0.0, [](double m, input_t x)->double{return m + x.realtime;}) / (double)cts2.size() : 0.0;
			double stepMeanC = ctsC.size() > 0 ? std::accumulate(ctsC.begin(), ctsC.end(), 0.0, [](double m, coinc_t x)->double{return m + x.realtime;}) / (double)ctsC.size() : 0.0;
			
			// Get dead time counts
			/*double time;
			if (cts1.size() != 0) {
				time = floor(cts1.front().realtime);
			} else if (cts2.size() != 0) {
				time = floor(cts2.front().realtime);
			} else {
				printf("No counts loaded! Continuing to next step!\n");
				continue;
			}*/

			//double corr1 = getDeadTimeSing( cts1, tempSta, tempEnd);
			//double corr2 = getDeadTimeSing( cts2, tempSta, tempEnd);
			double corrC = getDeadTimeCoinc(ctsC, tempSta, tempEnd);

			// Singles to coincidence conversion
			//double coin1 = -1.0; // -1.0 is the error value
			//double coin2 = -1.0;
			
			double pileup = puObj.CalculatePileup(ctsC);
			// TODO: Separate out corrC and Pileup corrections with a flag
	
			// Write out data!
			printf("\nData -----  Run: %d\nnumSteps: %d\ntempSta: %f\ntempEnd: %f\nstepMean1: %f\ncts1: %lu\ncorr1: %f\ncoin1: %f\nstepMean2: %f\ncts2: %lu\ncorr2: %f\ncoin2: %f\nstepMeanC: %f\nctsC: %lu\ncorrC: %f\n",
					runNo, numSteps, tempSta, tempEnd,
					stepMean1, cts1Size, corr1, coin1, stepMean2, cts2Size, corr2, coin2,
					stepMeanC, ctsCSize, corrC+pileup);
			fprintf(outfileSing,"%d,%d,%f,%f,%f,%lu,%f,%f,%f,%lu,%f,%f\n",
					runNo, numSteps, tempSta, tempEnd,
					stepMean1, cts1Size, corr1, coin1, stepMean2, cts2Size, corr2, coin2);
			fprintf(outfileCoin,"%d,%d,%f,%f,%f,%lu,%f\n",
					runNo, numSteps, tempSta, tempEnd, stepMeanC, ctsCSize, corrC+pileup);			
			//printf("PU = %f\n",pileup);
			pileupTot += pileup;
			// Increment start and end of step
			if (jj != 0) { // First slice has the movement:
				tempSta = tempEnd;
				tempEnd += sliceTime;
			} else {
				tempEnd = tempSta+sliceTime;
				tempSta = dagStop.at(sC);
				sC += 1;
				while(tempEnd < tempSta + 1.0) {
					tempEnd += sliceTime;
				}
			}				
		}
		numSteps += 1;		
	}
	fprintf(outfilePU,"%d,%f,%e,%f,",runNo,*dagSteps.begin() - (fillEnd+50),tau,puObj.GetAmplitude());
	fprintf(outfilePU,"%f\n",pileupTot);
	fclose(outfileSing);
	fclose(outfileCoin);	

}

//----------------------------------------------------------------------
void activeCleaner(Run* mcs1, Run* mcs2, EMS* ems0, const char* filename) {
	// Load tagbit data: FillEnd, H-GX, and Dagger Movement
	int runNo = mcs1->getRunNo();
	int cMode = mcs1->getCoincMode();
	double fillEnd = getFillEnd(mcs1, mcs2);
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2, fillEnd);
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
			
	// Make sure the run actually loaded...
	if(beamHits.empty() || (beamHits.end() <= beamHits.begin())) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", runNo);
		return;
	}
	if(dagSteps.empty() || (dagSteps.end() <= dagSteps.begin())) { 
		fprintf(stderr, "Skipping run %d for no Dagger Movement\n", runNo);
		return;
	}
	// Now find the background timing
	double countEnd;
	double bkgStart;
	if (dagSteps.size() > 2) { 
		countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
		bkgStart = mcs1->getTagBitEvt(1<<4, *(dagSteps.end()-2), 1) + 5.0;// 1<<4 is Cat Door, add extra 5 seconds to be safe
	} else if (dagSteps.size() > 1) { // For safety don't call if only one dagger step
		bkgStart = *(dagSteps.end()-1) - 40;
		fprintf(stderr, "\nUsing LAST 40s for background!!!!!!!\n");
	} else {
		fprintf(stderr,"Skipping run %d for non-existent dagger timings!\n", runNo);
		return;
	}
	double bkgEnd = *(dagSteps.end()-1);
	if (bkgStart < 5.0 || bkgStart >= bkgEnd) { // tag = -1.0 for errors (we've added 5s for stability.)
		fprintf(stderr, "Skipping run %d for background timing glitch\n", runNo);
		return;
	}
		
	// Lastly look at the cleaning time and the peak 1 time
	double pk1Start = *(dagSteps.begin());
	double pk1End;
	if (dagSteps.size() > 1) {
		pk1End   = *(dagSteps.begin()+1);
	} else {
		pk1End = pk1Start + 40;
		fprintf(stderr,"\nUsing 40 for peak 1 timing!!!!!!!\n");
	}
	double cleanEnd;
	
	double acUp    = mcs1->getTagBitEvt(1<<0, fillEnd + 1.0, 0);
	double gcUp    = mcs2->getTagBitEvt(1<<2, fillEnd + 1.0, 0);
	if ((acUp <= pk1Start) && (acUp >= 0.0)) { // Pick whichever one of these is 
		cleanEnd = acUp;
	} else if ((gcUp <= pk1Start) && (gcUp >= 0.0)) {
		cleanEnd = gcUp;
	} else {
		cleanEnd = fillEnd + 50.0;
		fprintf(stderr,"\nUsing 50s for cleaning time!!!!!!\n");
	}
	if (cleanEnd + 2.0 > pk1Start) { // Cleaner movement is glitchy in 2018
		cleanEnd = fillEnd + 50.0;
		fprintf(stderr,"\nUsing 50s for cleaning time!!!!!!\n");
	}
	
	// Make sure all the timings are in the right order:
	if ((cleanEnd < fillEnd) || (pk1Start < cleanEnd) || (pk1End < pk1Start) || (bkgStart < pk1End) || (bkgEnd < bkgStart)) {
		fprintf(stderr,"Skipping Run %d for inconisitent timings!", runNo);
		return;
	}
		
	// Load active cleaner coincidence structure
	std::vector<measurement> weightMon;
	std::vector<input_t> detCts;
	int detCh = (6 - 5) + mcs2->getMCSOff(); // hardcoding in just AC
	// Assuming active cleaner is ~2x the dagger window stuff
	double acCoinWindow = mcs1->getCoincWindow() * 2.0;
	double acPEWindow   = mcs1->getPeSumWindow() * 2.0;
	int acPESum = ceil(mcs1->getPeSum() / 2.0); 
	
	// Now just do vectors in the different time windows of interest
	std::vector<input_t> loadCts = mcs2->getCounts(
		[](input_t x)->input_t{return x;},
		[detCh, fillEnd](input_t x)->bool{return (x.ch == detCh && x.realtime < fillEnd);}
	);
	std::vector<coinc_t> loadCoincs = getSelfCoincs(loadCts,acCoinWindow,acPEWindow,acPESum);
	
	std::vector<input_t> cleanCts = mcs2->getCounts(
		[](input_t x)->input_t{return x;},
		[detCh, fillEnd, cleanEnd](input_t x)->bool{return (x.ch == detCh && x.realtime > fillEnd && x.realtime < cleanEnd);}
	);
	std::vector<coinc_t> cleanCoincs = getSelfCoincs(cleanCts,acCoinWindow,acPEWindow,acPESum);
	
	std::vector<input_t> holdCts = mcs2->getCounts(
		[](input_t x)->input_t{return x;},
		[detCh, cleanEnd, pk1Start](input_t x)->bool{return (x.ch == detCh && x.realtime > cleanEnd && x.realtime < pk1Start);}
	);
	std::vector<coinc_t> holdCoincs = getSelfCoincs(holdCts,acCoinWindow,acPEWindow,acPESum);
	
	std::vector<input_t> pk1Cts = mcs2->getCounts(
		[](input_t x)->input_t{return x;},
		[detCh, pk1Start, pk1End](input_t x)->bool{return (x.ch == detCh && x.realtime > pk1Start && x.realtime < pk1End);}
	);
	std::vector<coinc_t> pk1Coincs = getSelfCoincs(pk1Cts,acCoinWindow,acPEWindow,acPESum);
	
	std::vector<input_t> bkgCts = mcs2->getCounts(
		[](input_t x)->input_t{return x;},
		[detCh, bkgStart, bkgEnd](input_t x)->bool{return (x.ch == detCh && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	std::vector<coinc_t> bkgCoincs = getSelfCoincs(bkgCts,acCoinWindow,acPEWindow,acPESum);
	
	// And do output file
	FILE* outfileAC;
	outfileAC = fopen("acCoincs.csv","a");
	printf("Outputting active cleaner data to acCoincs.csv!\n");

	printf("\nActive Cleaner Data -----Run: %d\ntdTime: %f\nloadCts: %lu\ncleanTime: %f\ncleanCts:%lu\nholdTime:%f\nholdCts:%lu\npk1Time:%f\npk1Cts:%lu\nbkgTime:%f\nbkgCts:%lu\n",
			runNo, fillEnd, loadCoincs.size(),
			cleanEnd-fillEnd,cleanCoincs.size(),
			pk1Start-cleanEnd,holdCoincs.size(),
			pk1End-pk1Start,pk1Coincs.size(),
			bkgEnd-bkgStart,bkgCoincs.size());
	fprintf(outfileAC,"%d,%f,%lu,%f,%lu,%f,%lu,%f,%lu,%f,%lu\n",
			runNo, fillEnd, loadCoincs.size(),
			cleanEnd-fillEnd,cleanCoincs.size(),
			pk1Start-cleanEnd,holdCoincs.size(),
			pk1End-pk1Start,pk1Coincs.size(),
			bkgEnd-bkgStart,bkgCoincs.size());
	fclose(outfileAC);
}
//----------------------------------------------------------------------


//----------------------------------------------------------------------
// These functions might not really be useful anymore.
//----------------------------------------------------------------------
/* Normalized neutrons for Coincidence analysis */
void normNByDipCoinc(Run* mcs1, Run* mcs2, EMS* ems0, const char* filename) {
			
	// Load tagbit data: FillEnd, H-GX, and Dagger Movement
	int runNo = mcs1->getRunNo();
	int cMode = mcs1->getCoincMode();
	double fillEnd = getFillEnd(mcs1, mcs2);
	double tdTime = fillEnd;
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2, fillEnd);
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);
			
	// Make sure the run actually loaded...
	if(beamHits.empty() || (beamHits.end() <= beamHits.begin())) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", runNo);
		return;
	}
	if(dagSteps.empty() || (dagSteps.end() <= dagSteps.begin())) { 
		fprintf(stderr, "Skipping run %d for no Dagger Movement\n", runNo);
		return;
	}
	
	// load time constant data
	std::vector<int> rawDet;
	std::vector<double> kappa; 
	std::vector<double> kappaE;
	FILE* tc_file = fopen(filename, "r");
	if (tc_file != NULL) {
		char line[48];
		while (fgets(line,48,tc_file)) {
			// raw detectors coming from python script as 1-10. 
			const char* minRun = strtok(line,",");
			const char* maxRun = strtok(NULL,",");
			if ((std::stoi(minRun) <= runNo) && (std::stoi(maxRun) > runNo)) {
				const char* detF = strtok(NULL,",");
				rawDet.push_back(std::stoi(detF));
				kappa.push_back(std::stod(strtok(NULL,",")));
				kappaE.push_back(std::stod(strtok(NULL,",")));
			}
		}
		fclose(tc_file);
	}
	if (rawDet.size() == 0) { // If there's a problem in loading tc_file assume {N, 0.0,0.0}
		//printf("No TC file found, assuming 70 s fit!\n");
		for (int i = 1; i < 11; i++) {
			rawDet.push_back(i);
			kappa.push_back(0.0); // No tc_file means no weighting
			kappaE.push_back(0.0);
		}
	}	
	// Load a vector of weighted monitors.
	Run* dagMCS;
	std::vector<measurement> weightMon;
	std::vector<input_t> detCts;
	for (int ii = 0; ii < rawDet.size(); ii++) {
		int detCh;
		// need to convert detectors 1-10 to the right MCS box
		if (ii < 5) {
			dagMCS = mcs1;
			detCh = rawDet[ii] + dagMCS->getMCSOff();
		} else {
			dagMCS = mcs2;
			detCh = (rawDet[ii]-5) + dagMCS->getMCSOff();
		}
		detCts = dagMCS->getCounts(
			[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
			[detCh, fillEnd](input_t x)->bool{return x.ch == detCh && x.realtime < fillEnd;}
		);
	
		weightMon.push_back(expWeightMonVect(detCts, kappa[ii]));
	}
	
	// Determine the threshold
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	// Find background times
	// During production running, the last 50s should have TD open, dagger down bkgs
	double bkgStart;
	double bkgEnd = *(dagSteps.end()-1);
	if (dagSteps.size() > 2) { 
		bkgStart = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
	} else {
		bkgStart = *(dagSteps.end()-1) - 50;
		fprintf(stderr, "\nUsing LAST 50s for background!!!!!!!\n");
	}
	
	if (bkgStart < 0.0 || bkgStart >= bkgEnd) {
		fprintf(stderr, "Skipping run %d for background timing glitch\n", runNo);
		return;
	}
	
	// Create vector of background coincidence and singles counts
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	std::vector<input_t> bkgDag1;
	std::vector<input_t> bkgDag2;
	std::vector<coinc_t> bkgDagC;
	bkgDag1 = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c1 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	bkgDag2 = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c2 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	bkgDagC = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[bkgStart, bkgEnd](coinc_t x)->bool{return (x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	
	// Output 2 .csv files that we will use in NormAnalyzer-Coinc.py
	// Start by outputting monitor and backgrounds
	FILE* outfileDet;
	outfileDet = fopen("normNByDipCoinc-det.csv", "a");
	printf("Outputting monitor/background data to normNByDipCoinc-det.csv!\n");
	printf("\nMonitor Data -----Run: %d\ntdTime: %f\nholdEnd: %f\bkgDag1: %lu\nbkgDag2: %lu\nbkgDagC: %lu\nbkgStart: %f\nbkgEnd: %f\n",
			runNo, tdTime, *dagSteps.begin(), bkgDag1.size(),bkgDag2.size(),bkgDagC.size(), bkgStart, bkgEnd);
	fprintf(outfileDet,"%d,%f,%f,%lu,%lu,%lu,%f,%f",
			runNo, tdTime, *dagSteps.begin(), bkgDag1.size(),bkgDag2.size(),bkgDagC.size(), bkgStart, bkgEnd);			
	for (int jt = 0; jt < weightMon.size(); jt++ ) {
		printf("Monitor %d: %f (%f)\n",jt,weightMon[jt].val,weightMon[jt].err);
		fprintf(outfileDet,",%f,%f",weightMon[jt].val,weightMon[jt].err);
	}
	fprintf(outfileDet,"\n");
	fclose(outfileDet);
	
	// Now output dagger detector data
	FILE* outfileDag;
	outfileDag = fopen("normNByDipCoinc-dag.csv","a");
	
	// Separating dips into 10s slices (arb.)
	int numSteps = 0;
	double sliceTime = 10.0;
		
	// Loop through each dagger step and find when photon events hit
	for(auto stepIt = dagSteps.begin()+1; stepIt < dagSteps.end(); stepIt++) {
		
		// Time for this step
		double startTime = *(stepIt-1);
		double endTime = *(stepIt);
						
		// Cut off the background period in the last 50s.
		if(stepIt == dagSteps.end()-1) { 
			endTime = bkgStart;
		}
		
		// Divide each step into 10 second slices
		double tempSta = startTime;
		double tempEnd = startTime + sliceTime;
		for (int jj = 0; jj < (int)((endTime - startTime) / sliceTime) + 1; jj++) {
			
			// Make sure we stop at the right step
			if (tempEnd > endTime) {
				tempEnd = endTime;
			}
			// 0.2s timing slop
			if ((tempSta + 0.2) > tempEnd) { 
				break;
			}
			
			// Read in data (coincidence)
			std::vector<coinc_t> ctsC = dagMCS->getCoincCounts(
				[](coinc_t x)->coinc_t{return x;},
				[tempSta, tempEnd](coinc_t x)->bool{return (x.realtime > tempSta && x.realtime < tempEnd);}
			);
			// Find Mean Arrival Time
			double stepMeanC = ctsC.size() > 0 ?
				std::accumulate(ctsC.begin(), ctsC.end(), 0.0, [](double m, coinc_t x)->double{return m + x.realtime;}) / (double)ctsC.size()
				: 0.0;
			
			// Get dead time counts (shouldn't matter since we're choosing ideal rate-dep values)
			double corrC = getDeadTimeCoinc(ctsC, tempSta, tempEnd);

			printf("\nData -----  Run: %d\numSteps: %d\ntempSta: %f\ntempEnd: %f\nstepMeanC: %f\nctsC: %lu\ncorrC: %f\n",
					runNo, numSteps, tempSta, tempEnd, stepMeanC, ctsC.size(), corrC);
			fprintf(outfileDag,"%d,%d,%f,%f,%f,%lu,%f\n",
					runNo, numSteps, tempSta, tempEnd, stepMeanC, ctsC.size(), corrC);
	
			// Increment start and end
			tempSta += sliceTime;
			tempEnd += sliceTime;
		}
		numSteps += 1;
	}
	fclose(outfileDag);
}


void extractObservables(Run* mcs1, Run* mcs2, const char* filename) {
	// Begin by loading things
	int cMode = mcs1->getCoincMode();
	int runNo = mcs1->getRunNo();
	double fillEnd = getFillEnd(mcs1, mcs2);
	double tdTime = fillEnd;
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2, fillEnd);
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);
	
	// Make sure the run actually loaded.
	if(beamHits.empty() || (beamHits.end() <= beamHits.begin())) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", mcs1->getRunNo());
		return;
	}
	if(beamHits.empty() || (beamHits.end() <= beamHits.begin())) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", mcs1->getRunNo());
		return;
	}
	
	// If there are no dagger movements, skip the run (it was probably killed due to bad fill) 
	if (dagSteps.end() == dagSteps.begin()) { 
		return;
	}
	
	// load time constant data
	FILE* tc_file = fopen(filename, "r");
	char line[48];
	std::vector<int> rawDet;
	std::vector<double> kappa; 
	std::vector<double> kappaE;
	while (fgets(line,48,tc_file)) {
		// raw detectors coming from python script as 1-10. 
		const char* minRun = strtok(line,",");
		const char* maxRun = strtok(NULL,",");
		if ((std::stoi(minRun) <= runNo) && (std::stoi(maxRun) > runNo)) {
			const char* detF = strtok(NULL,",");
			rawDet.push_back(std::stoi(detF));
			kappa.push_back(std::stod(strtok(NULL,",")));
			kappaE.push_back(std::stod(strtok(NULL,",")));
		}
	}
	fclose(tc_file);
		
	// Load a vector of weighted monitors -- only works for 2017 
	std::vector<input_t> oldCts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 3 && x.realtime < fillEnd;}
	);
	std::vector<input_t> stpCts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 5 && x.realtime < fillEnd;}
	);
	std::vector<input_t> barCts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 4 && x.realtime < fillEnd;}
	);
	
	// Raw weighting for the monitors
	int oldMon = (int)oldCts.size();
	int barMon = (int)barCts.size();
	int stpMon = (int)stpCts.size();
		
	// Calculate exponential weighting for these three monitors
	measurement weightBar = expWeightMonVect(barCts,kappa[4]);
	measurement weightOld = expWeightMonVect(oldCts,kappa[3]);
	measurement weightStp = expWeightMonVect(stpCts,kappa[5]);
	
	// Find background times
	// During production running, the last 50s should have TD open, dagger down bkgs
	double bkgStart;
	double bkgEnd = *(dagSteps.end()-1);
	
	// Find where the end of run happens and EOR backgrounds start.
	if (dagSteps.size() > 2) { 
		bkgStart = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
	} else {
		bkgStart = *(dagSteps.end()-1) - 50;
		fprintf(stderr, "\nUsing LAST 50s for background!!!!!!!\n");
	}
	
	if (bkgStart < 0.0 || bkgStart >= bkgEnd) {
		fprintf(stderr, "Skipping run %d for background timing glitch\n", runNo);
		return;
	}
	// Load counts on background detectors
	std::vector<input_t> bkgD1Cts;
	std::vector<input_t> bkgD2Cts;
	// Check high or low threshold
	if ((cMode == 1) || (cMode == 2)) {
		bkgD1Cts = mcs1->getCounts(
			[](input_t x)->input_t{return x;},
			[bkgStart,bkgEnd](input_t x)->bool{return (x.ch == 1 && x.realtime > bkgStart && x.realtime < bkgEnd);}
		);
		bkgD2Cts = mcs1->getCounts(
			[](input_t x)->input_t{return x;},
			[bkgStart,bkgEnd](input_t x)->bool{return (x.ch == 2 && x.realtime > bkgStart && x.realtime < bkgEnd);}
		);
	} else if ((cMode == 3) || (cMode == 4)) {		
		bkgD1Cts = mcs2->getCounts(
			[](input_t x)->input_t{return x;},
			[bkgStart,bkgEnd](input_t x)->bool{return (x.ch == 14 && x.realtime > bkgStart && x.realtime < bkgEnd);}
		);
		bkgD2Cts = mcs2->getCounts(
			[](input_t x)->input_t{return x;},
			[bkgStart,bkgEnd](input_t x)->bool{return (x.ch == 15 && x.realtime > bkgStart && x.realtime < bkgEnd);}
		);
	}
	std::vector<coinc_t> bkgCCts = mcs1->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[bkgStart, bkgEnd](coinc_t x)->bool{return (x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	
	// Convert vector counts to rates for bkg ratios
	int bkgD1 = (int)bkgD1Cts.size();
	int bkgD2 = (int)bkgD2Cts.size();
	int bkgDC = (int)bkgCCts.size();
	
	// Output a .csv file that we can use
	FILE* outfile;
	outfile = fopen("ExtractedObservables.csv", "a");
	printf("Outputting data to ExtractedObservables.csv!\n");
	
	int ii = 0;
	int dipCts1 [3];
	int dipCts2 [3];
	int dipCtsC [3];
	// Loop through each dagger step and find when photon events hit
	//for(auto stepIt = dagSteps.begin()+1; stepIt < dagSteps.end(); stepIt++) {
	for(auto stepIt = dagSteps.begin()+1; stepIt < dagSteps.end(); stepIt++) {
		if (ii > 2) {break;}
		double startTime = *(stepIt-1);
		double endTime = *(stepIt);
					
		// Cut off the background period.
		if(stepIt == dagSteps.end()-1) { 
			endTime = bkgStart;
		}
		if (startTime >= endTime) {
			printf("Stopping due to timing bug!\n");
			return;
		}
					
		// Load counts, and buffer
		std::vector<input_t> d1Cts;
		std::vector<input_t> d2Cts;
		std::vector<coinc_t> dCCts;	
		
		if ((cMode == 1) || (cMode == 2)) {
			d1Cts = mcs1->getCounts(
				[](input_t x)->input_t{return x;},
				[startTime, endTime](input_t x)->bool{return (x.ch == 1 && x.realtime > startTime && x.realtime < endTime);}
			);
			d2Cts = mcs1->getCounts(
				[](input_t x)->input_t{return x;},
				[startTime, endTime](input_t x)->bool{return (x.ch == 2 && x.realtime > startTime && x.realtime < endTime);}
			);
			
			dCCts = mcs1->getCoincCounts(
				[](coinc_t x)->coinc_t{return x;},
				[startTime, endTime](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
			);
			
		// High threshold
		} else if ((cMode == 3) || (cMode == 4)) {
			d1Cts = mcs2->getCounts(
				[](input_t x)->input_t{return x;},
				[startTime, endTime](input_t x)->bool{return (x.ch == 14 && x.realtime > startTime && x.realtime < endTime);}
			);
			d2Cts = mcs2->getCounts(
				[](input_t x)->input_t{return x;},
				[startTime, endTime](input_t x)->bool{return (x.ch == 15 && x.realtime > startTime && x.realtime < endTime);}
			);
			dCCts = mcs2->getCoincCounts(
				[](coinc_t x)->coinc_t{return x;},
				[startTime, endTime](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
			);
			
		}	
		dipCts1[ii] = (int)d1Cts.size();
		dipCts2[ii] = (int)d2Cts.size();
		dipCtsC[ii] = (int)dCCts.size();
		ii+=1;
			
	}

	fprintf(outfile,"%d,%d,%d,%d,%d,%f,%f,%f,%f,%d,%d,%d,%d,%f,%f,%f,%f,%d,%d,%d,%d,%f,%f,%f,%f,%d,%d,%d,%f,%f,%f\n",
			mcs1->getRunNo(), 
			dipCtsC[0], dipCtsC[1],dipCtsC[2],bkgDC,
			sqrt((float)dipCtsC[0]),sqrt((float)dipCtsC[1]),sqrt((float)dipCtsC[2]),sqrt(bkgDC),
			dipCts1[0], dipCts1[1],dipCts1[2],bkgD1,
			sqrt((float)dipCts1[0]),sqrt((float)dipCts1[1]),sqrt((float)dipCts1[2]),sqrt(bkgD1),
			dipCts2[0], dipCts2[1],dipCts2[2],bkgD2,
			sqrt((float)dipCts2[0]),sqrt((float)dipCts2[1]),sqrt((float)dipCts2[2]),sqrt(bkgD2),
			oldMon,stpMon,barMon,weightOld.val,weightStp.val,weightBar.val);

	fclose(outfile);
	
}


/* Generate normalized fill data only (does not re-run coincidence or 
 * singles detector unload counts) */
void normNByDipDet(Run* mcs1, Run* mcs2, EMS* ems0, const char* filename) {
	
	// Load tagbit data: FillEnd, H-GX, and Dagger Movement
	int runNo = mcs1->getRunNo();
	int cMode = mcs1->getCoincMode();
	double fillEnd = getFillEnd(mcs1, mcs2);
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2, fillEnd);
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
			
	// Make sure the run actually loaded...
	if(beamHits.empty() || (beamHits.end() <= beamHits.begin())) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", runNo);
		return;
	}
	if(dagSteps.empty() || (dagSteps.end() <= dagSteps.begin())) { 
		fprintf(stderr, "Skipping run %d for no Dagger Movement\n", runNo);
		return;
	}

	// load time constant data
	std::vector<int> rawDet;
	std::vector<double> kappa; 
	std::vector<double> kappaE;
	FILE* tc_file = fopen(filename, "r");
	if (tc_file != NULL) {
		char line[48];
		while (fgets(line,48,tc_file)) {
			// raw detectors coming from python script as 1-10. 
			const char* minRun = strtok(line,",");
			const char* maxRun = strtok(NULL,",");
			if ((std::stoi(minRun) <= runNo) && (std::stoi(maxRun) > runNo)) {
				const char* detF = strtok(NULL,",");
				rawDet.push_back(std::stoi(detF));
				kappa.push_back(std::stod(strtok(NULL,",")));
				kappaE.push_back(std::stod(strtok(NULL,",")));
			}
		}
		fclose(tc_file);
	}
	std::vector<double> mBkg;
	std::vector<double> mBkgE;
	const char* monBkgLoc  = std::getenv("MONBKG_LOC"); 
	FILE* mb_file = fopen(monBkgLoc,"r");
	if (mb_file != NULL) {
		char line[48];
		while (fgets(line,48,mb_file)) {
			// Monitor Backgrounds coming from python script as 1-10.
			const char* minRun = strtok(line,",");
			const char* maxRun = strtok(NULL,",");
			if ((std::stoi(minRun) <= runNo) && (std::stoi(maxRun) > runNo)) {
				const char* detF = strtok(NULL,",");
				mBkg.push_back(std::stod(strtok(NULL,",")));
				mBkgE.push_back(std::stod(strtok(NULL,",")));
			}
		}
		fclose(mb_file);
	}

	if (rawDet.size() == 0) { // If there's a problem in loading tc_file assume {N, 0.0,0.0}
		//printf("No TC file found, assuming 70 s fit!\n");
		for (int i = 1; i < 11; i++) {
			rawDet.push_back(i);
			kappa.push_back(0.0); // No tc_file means no weighting
			kappaE.push_back(0.0);
		}
	}
	if (rawDet.size() != mBkg.size()) {
		mBkg.clear();
		mBkgE.clear();
		for (size_t i = 0; i < rawDet.size(); i++){
			mBkg.push_back(0.0);
			mBkgE.push_back(0.0);
		}
	}
		
	// Find background times
	// During production running, the last 50s should have TD open, dagger down bkgs
	/*double bkgStart;
	double bkgEnd = *(dagSteps.end()-1);
	if (dagSteps.size() > 2) { 
		bkgStart = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
	} else {
		bkgStart = *(dagSteps.end()-1) - 50;
		fprintf(stderr, "\nUsing LAST 50s for background!!!!!!!\n");
	}
	
	if (bkgStart < 0.0 || bkgStart >= bkgEnd) {
		fprintf(stderr, "Skipping run %d for background timing glitch\n", runNo);
		return;
	}*/
	double bkgStart;
	double bkgEnd = *(dagSteps.end()-1);
	if (dagSteps.size() > 2) { 
	//	lastDip = *(dagSteps.end() - 2);
	//	countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
		bkgStart = mcs1->getTagBitEvt(1<<4, *(dagSteps.end()-2), 1) + 5.0;// 1<<4 is Cat Door, add extra 5 seconds to be safe
	} else {
		bkgStart = *(dagSteps.end()-1) - 40;
		fprintf(stderr, "\nUsing LAST 40s for background!!!!!!!\n");
	}
	
	if (bkgStart < 5.0 || bkgStart >= bkgEnd) { // tag = -1.0 for errors (we've added 5s for stability.)
		fprintf(stderr, "Skipping run %d for background timing glitch\n", runNo);
		return;
	}
	
	
	double bkgMStart = *(dagSteps.end()-1) - 50; // Hardcoding last 50s.
	
	// Load a vector of weighted monitors.
	Run* dagMCS;
	// Determine the threshold -- not really needed here
	if ((mcs1->getCoincMode() % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
		
	//std::vector<long> weightMon; // ignore uncertainties for now
	std::vector<measurement> weightMon;
	std::vector<long> bkgMon;  //Bakground of monitor counts
	for (int ii = 0; ii < rawDet.size(); ii++) {
		int detCh;
		// need to convert detectors 1-10 to the right MCS box
		if (ii < 5) {
			dagMCS = mcs1;
			detCh = rawDet[ii] + dagMCS->getMCSOff();
		} else {
			dagMCS = mcs2;
			detCh = (rawDet[ii]-5) + dagMCS->getMCSOff();
		}
		std::vector<input_t> detCts = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			//[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
			[detCh, fillEnd](input_t x)->bool{return x.ch == detCh && x.realtime < fillEnd;}
		);
		//printf("%lu,%lu\n",detCts,
		//printf("%lu,%lu\n",detCts,
		std::vector<input_t> bkgCts = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[detCh, bkgMStart,bkgEnd](input_t x)->bool{return x.ch == detCh && (bkgMStart <= x.realtime && x.realtime < bkgEnd);}
		);
		//if (rawDet[ii] != 4 && runNo >= 10988) { // Putting in self-integration for detector 4
		weightMon.push_back(expWeightMonVectBkg(detCts, kappa[ii],16*NANOSECOND,(double)bkgCts.size()/(bkgEnd-bkgMStart),fillEnd));
		//	weightMon.push_back((double)detCts.size()); // Raw monitor counts
			bkgMon.push_back(bkgCts.size());
		/*} else {
			std::vector<coinc_t> parseCts = getSelfCoincs(detCts,1500,0,2);
			std::vector<coinc_t> parseBkg = getSelfCoincs(bkgCts,1500,0,2);
			weightMon.push_back(parseCts.size());
			bkgMon.push_back(parseBkg.size());
			//bkgMon.push_back(bkgCts.size());
		}*/	

	}
	
	/*FILE* outfileDet;
	outfileDet = fopen("detectorAndBkg.csv","a");
	printf("Outputting monitor+background data to detectorAndBkg.csv!\n");
	printf("\nMonitor Data -- Run: %d\nHold: %f\nBkgLen: %f\n", 
			runNo, *dagSteps.begin() - (fillEnd+50), bkgEnd - bkgStart);
	fprintf(outfileDet,"%d,%f,%f",runNo,*dagSteps.begin()-(fillEnd+50),bkgEnd-bkgStart);
	for (size_t jt = 0; jt < weightMon.size(); jt++ ) {
		printf("Monitor %d: %u (%u)\n",jt,weightMon[jt],bkgMon[jt]);
		fprintf(outfileDet,",%u,%u",weightMon[jt],bkgMon[jt]);
	}
	fprintf(outfileDet,"\n");*/
	std::vector<input_t> bkgDag1 = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c1 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	std::vector<input_t> bkgDag2 = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c2 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	std::vector<coinc_t> bkgDagC = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[bkgStart, bkgEnd](coinc_t x)->bool{return (x.realtime > bkgStart && x.realtime < bkgEnd);}
	);		
	
		// Output 2 .csv files that we will use in NormAnalyzer-Coinc.py
	// Start by outputting monitor and backgrounds
	FILE* outfileDet;
	outfileDet = fopen("normNByDip-det_B.csv", "a");
	printf("Outputting monitor/background data to normNByDip-det_B.csv!\n");
	//printf("\nMonitor Data -----Run: %d\ntdTime: %f\nholdEnd: %f\nbkgDag1: %lu\nbkgDag2: %lu\nbkgDagC: %lu\nbkgStart: %f\nbkgEnd: %f\n",
	//		runNo, fillEnd, *dagSteps.begin(), bkgDag1Size,bkgDag2Size,bkgDagCSize, bkgStart, bkgEnd);
	fprintf(outfileDet,"%d,%f,%f,%lu,%lu,%lu,%f,%f",
			runNo, fillEnd, *dagSteps.begin(), bkgDag1.size(),bkgDag2.size(),bkgDagC.size(), bkgStart, bkgEnd);			
	for (int jt = 0; jt < weightMon.size(); jt++ ) {
		printf("Monitor %d: %f (%f)\n",jt,weightMon[jt].val,weightMon[jt].err);
		fprintf(outfileDet,",%f,%f",weightMon[jt].val,weightMon[jt].err);
	}
	fprintf(outfileDet,"\n");
	fclose(outfileDet);
	
	
	// Load background counts
	/*std::vector<input_t> bkgDag1;
	std::vector<input_t> bkgDag2;
	std::vector<coinc_t> bkgDagC;
	bkgDag1 = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c1 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	bkgDag2 = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c2 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	bkgDagC = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[bkgStart, bkgEnd](coinc_t x)->bool{return (x.realtime > bkgStart && x.realtime < bkgEnd);}
	);*/			
			
	// Output 2 .csv files that we will use in NormAnalyzer-Coinc.py
	// Start by outputting monitor and backgrounds
	/*FILE* outfileDet;
	outfileDet = fopen("normNByDip-det.csv", "a");
	printf("Outputting monitor/background data to normNByDip-det.csv!\n");
	printf("\nMonitor Data -----Run: %d\ntdTime: %f\nholdEnd: %f\nbkgDag1: %lu\nbkgDag2: %lu\nbkgDagC: %lu\nbkgStart: %f\nbkgEnd: %f\n",
			runNo, fillEnd, *dagSteps.begin(), bkgDag1.size(),bkgDag2.size(),bkgDagC.size(), bkgStart, bkgEnd);
	fprintf(outfileDet,"%d,%f,%f,%lu,%lu,%lu,%f,%f",
			runNo, fillEnd, *dagSteps.begin(), bkgDag1.size(),bkgDag2.size(),bkgDagC.size(), bkgStart, bkgEnd);			
	for (int jt = 0; jt < weightMon.size(); jt++ ) {
		printf("Monitor %d: %f (%f)\n",jt,weightMon[jt].val,weightMon[jt].err);
		fprintf(outfileDet,",%f,%f",weightMon[jt].val,weightMon[jt].err);
	}
	fprintf(outfileDet,"\n");
	fclose(outfileDet);*/
	
}



/* Normalized neutrons for Singles analysis -- this will do arbitrary 
 * photon thresholds.*/
void normNByDipSing(Run* mcs1, Run* mcs2, EMS* ems0, const char* filename) {
	
	// Load tagbit data: FillEnd, H-GX, and Dagger Movement
	int runNo = mcs1->getRunNo();
	int cMode = mcs1->getCoincMode();
	double fillEnd = getFillEnd(mcs1, mcs2);
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2, fillEnd);
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
			
	// Make sure the run actually loaded...
	if(beamHits.empty() || (beamHits.end() <= beamHits.begin())) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", runNo);
		return;
	}
	if(dagSteps.empty() || (dagSteps.end() <= dagSteps.begin())) { 
		fprintf(stderr, "Skipping run %d for no Dagger Movement\n", runNo);
		return;
	}

	// load time constant data
	std::vector<int> rawDet;
	std::vector<double> kappa; 
	std::vector<double> kappaE;
	FILE* tc_file = fopen(filename, "r");
	if (tc_file != NULL) {
		char line[48];
		while (fgets(line,48,tc_file)) {
			// raw detectors coming from python script as 1-10. 
			const char* minRun = strtok(line,",");
			const char* maxRun = strtok(NULL,",");
			if ((std::stoi(minRun) <= runNo) && (std::stoi(maxRun) > runNo)) {
				const char* detF = strtok(NULL,",");
				rawDet.push_back(std::stoi(detF));
				kappa.push_back(std::stod(strtok(NULL,",")));
				kappaE.push_back(std::stod(strtok(NULL,",")));
			}
		}
		fclose(tc_file);
	}
	if (rawDet.size() == 0) { // If there's a problem in loading tc_file assume {N, 0.0,0.0}
		//printf("No TC file found, assuming 70 s fit!\n");
		for (int i = 1; i < 11; i++) {
			rawDet.push_back(i);
			kappa.push_back(0.0); // No tc_file means no weighting
			kappaE.push_back(0.0);
		}
	}
			
	// Load a vector of weighted monitors.
	Run* dagMCS;
	std::vector<measurement> weightMon;
	std::vector<input_t> detCts;
	for (int ii = 0; ii < rawDet.size(); ii++) {
		int detCh;
		// need to convert detectors 1-10 to the right MCS box
		if (ii < 5) {
			dagMCS = mcs1;
			detCh = rawDet[ii] + dagMCS->getMCSOff();
		} else {
			dagMCS = mcs2;
			detCh = (rawDet[ii]-5) + dagMCS->getMCSOff();
		}
		detCts = dagMCS->getCounts(
			[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
			[detCh, fillEnd](input_t x)->bool{return x.ch == detCh && x.realtime < fillEnd;}
		);
	
		weightMon.push_back(expWeightMonVectDT(detCts, kappa[ii],16*NANOSECOND,fillEnd));
	}
	
	// Find background times
	// During production running, the last 50s should have TD open, dagger down bkgs
	double countEnd;
	double bkgStart;
	double bkgEnd = *(dagSteps.end()-1);
	if (dagSteps.size() > 2) { 
		countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
		bkgStart = mcs1->getTagBitEvt(1<<4, *(dagSteps.end()-2), 1) + 5.0;// 1<<4 is Cat Door, add extra 5 seconds to be safe
	} else {
		bkgStart = *(dagSteps.end()-1) - 40;
		fprintf(stderr, "\nUsing LAST 40s for background!!!!!!!\n");
	}
	
	if (bkgStart < 5.0 || bkgStart >= bkgEnd) { // tag = -1.0 for errors (we've added 5s for stability.)
		fprintf(stderr, "Skipping run %d for background timing glitch\n", runNo);
		return;
	}
	
	// Here we're going to do 
	// Determine which dagger pair we want to use
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	// Load background counts
	std::vector<input_t> bkgTmp1;
	std::vector<input_t> bkgTmp2;
	std::vector<coinc_t> bkgDagC;
	bkgTmp1 = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c1 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	bkgTmp2 = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c2 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	// Assuming active cleaner is ~2x the dagger window stuff
	double coinWindow = mcs1->getCoincWindow(); // Initial window
	double PEWindow   = mcs1->getPeSumWindow(); // Telescope window
	int PESum = ceil(mcs1->getPeSum() / 2.0);  // Half the photon thr.
	
	// Need to self-coincidence these two counts
	std::vector<coinc_t> bkgDag1 = getSelfCoincs(bkgTmp1, coinWindow, PEWindow, PESum);
	std::vector<coinc_t> bkgDag2 = getSelfCoincs(bkgTmp2, coinWindow, PEWindow, PESum);
	
	bkgDagC = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[bkgStart, bkgEnd](coinc_t x)->bool{return (x.realtime > bkgStart && x.realtime < bkgEnd);}
	); // Should use this on 5 or 6 so this is cross PMT
			
	// Output 2 .csv files that we will use in LifetimeAnalyzer python.
	// Start by outputting monitor and backgrounds
	FILE* outfileDet;
	outfileDet = fopen("normNByDip-det.csv", "a");
	printf("Outputting monitor/background data to normNByDip-det.csv!\n");
	printf("\nMonitor Data -----Run: %d\ntdTime: %f\nholdEnd: %f\nbkgDag1: %lu\nbkgDag2: %lu\nbkgDagC: %lu\nbkgStart: %f\nbkgEnd: %f\n",
			runNo, fillEnd, *dagSteps.begin(), bkgDag1.size(),bkgDag2.size(),bkgDagC.size(), bkgStart, bkgEnd);
	fprintf(outfileDet,"%d,%f,%f,%lu,%lu,%lu,%f,%f",
			runNo, fillEnd, *dagSteps.begin(), bkgDag1.size(),bkgDag2.size(),bkgDagC.size(), bkgStart, bkgEnd);			
	for (int jt = 0; jt < weightMon.size(); jt++ ) {
		printf("Monitor %d: %f (%f)\n",jt,weightMon[jt].val,weightMon[jt].err);
		fprintf(outfileDet,",%f,%f",weightMon[jt].val,weightMon[jt].err);
	}
	fprintf(outfileDet,"\n");
	fclose(outfileDet);
	
	// Now output dagger detector data
	FILE* outfileSing;
	outfileSing = fopen("normNByDipSing-dag.csv","a");
	FILE* outfileCoin;
	outfileCoin = fopen("normNByDipCoinc-dag.csv","a");
	
	// Separating dips into 10s slices (arb.)
	int numSteps = 0;
	double sliceTime = 10.0;
		
	// Loop through each dagger step and find when photon events hit
	for(auto stepIt = dagSteps.begin()+1; stepIt < dagSteps.end(); stepIt++) {
		
		// Time for this step
		double startTime = *(stepIt-1);
		double endTime = *(stepIt);
						
		// Cut off the background period in the last 50s.
		if(stepIt == dagSteps.end()-1) { 
			endTime = countEnd;
		}
		
		if (startTime >= endTime) {
			printf("Stopping due to timing bug!\n");
			return;
		}
		
		// Divide each step into 10 second slices
		double tempSta = startTime;
		double tempEnd = startTime + sliceTime;
		for (int jj = 0; jj < (int)((endTime - startTime) / sliceTime) + 1; jj++) {
			
			// Make sure we stop at the right step
			if (tempEnd > endTime) {
				tempEnd = endTime;
			}
			// 0.2s timing slop
			if ((tempSta + 0.2) > tempEnd) { 
				break;
			}
			
			// Load counts (singles + coincidence)
			std::vector<input_t> tmp1 = dagMCS->getCounts(
				[](input_t x)->input_t{return x;},
				[c1,tempSta,tempEnd](input_t x)->bool{return (x.ch == c1 && x.realtime > tempSta && x.realtime < tempEnd);}
			);
			std::vector<input_t> tmp2 = dagMCS->getCounts(
				[](input_t x)->input_t{return x;},
				[c2,tempSta,tempEnd](input_t x)->bool{return (x.ch == c2 && x.realtime > tempSta && x.realtime < tempEnd);}
			);
			std::vector<coinc_t> cts1 = getSelfCoincs(tmp1, coinWindow, PEWindow, PESum);
			std::vector<coinc_t> cts2 = getSelfCoincs(tmp2, coinWindow, PEWindow, PESum);
			
			std::vector<coinc_t> ctsC = dagMCS->getCoincCounts(
				[](coinc_t x)->coinc_t{return x;},
				[tempSta, tempEnd](coinc_t x)->bool{return (x.realtime > tempSta && x.realtime < tempEnd);}
			);
			
			// Find Mean Arrival Time
			double stepMean1 = cts1.size() > 0 ? std::accumulate(cts1.begin(), cts1.end(), 0.0, [](double m, coinc_t x)->double{return m + x.realtime;}) / (double)cts1.size() : 0.0;
			double stepMean2 = cts2.size() > 0 ? std::accumulate(cts2.begin(), cts2.end(), 0.0, [](double m, coinc_t x)->double{return m + x.realtime;}) / (double)cts2.size() : 0.0;
			double stepMeanC = ctsC.size() > 0 ? std::accumulate(ctsC.begin(), ctsC.end(), 0.0, [](double m, coinc_t x)->double{return m + x.realtime;}) / (double)ctsC.size() : 0.0;
			
			// Get dead time counts
			double time;
			if (cts1.size() != 0) {
				time = floor(cts1.front().realtime);
			} else if (cts2.size() != 0) {
				time = floor(cts2.front().realtime);
			} else {
				printf("No counts loaded! Continuing to next step!\n");
				continue;
			}

			// Singles hardcoded in s with instantaneous rates every s.
			double corr1 = getDeadTimeCoinc(cts1,tempSta,tempEnd);
			double corr2 = getDeadTimeCoinc(cts2,tempSta,tempEnd);
			double corrC = getDeadTimeCoinc(ctsC, tempSta, tempEnd);

			// Singles to coincidence conversion
			double coin1 = -1.0; // -1.0 is the error value
			double coin2 = -1.0;
			if (ctsC.size() > 0) { // divide by zero
				coin1 = 0.0;
				coin2 = 0.0;
				for (auto cIt = ctsC.begin(); cIt < ctsC.end(); cIt++) {
					coin1 += (double)cIt->pmt1;
					coin2 += (double)cIt->pmt2;
				}
				coin1 = coin1 / (double)ctsC.size();
				coin2 = coin2 / (double)ctsC.size();
			}
				
			// Write out data!
			printf("\nData -----  Run: %d\nnumSteps: %d\ntempSta: %f\ntempEnd: %f\nstepMean1: %f\ncts1: %lu\ncorr1: %f\ncoin1: %f\nstepMean2: %f\ncts2: %lu\ncorr2: %f\ncoin2: %f\nstepMeanC: %f\nctsC: %lu\ncorrC: %f\n",
					runNo, numSteps, tempSta, tempEnd,
					stepMean1, cts1.size(), corr1, coin1, stepMean2, cts2.size(), corr2, coin2,
					stepMeanC, ctsC.size(), corrC);

			fprintf(outfileSing,"%d,%d,%f,%f,%f,%lu,%f,%f,%f,%lu,%f,%f\n",
					runNo, numSteps, tempSta, tempEnd,
					stepMean1, cts1.size(), corr1, coin1, stepMean2, cts2.size(), corr2, coin2);
			fprintf(outfileCoin,"%d,%d,%f,%f,%f,%lu,%f\n",
					runNo, numSteps, tempSta, tempEnd, stepMeanC, ctsC.size(), corrC);			

			// Increment start and end of step
			tempSta += sliceTime;
			tempEnd += sliceTime;
		}
		numSteps += 1;		
	}
	fclose(outfileSing);
	fclose(outfileCoin);
}
