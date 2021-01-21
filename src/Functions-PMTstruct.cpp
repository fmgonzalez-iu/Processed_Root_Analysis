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
	Author: Frank M. Gonzalez

	These functions look at the structure of events hitting the PMTs.
	A lot of these are outputting timing ROOT structures.
	
	bkgTTNE and getTtne are deprecated time to next event histogram generators.
	
------------------------------------------------------------------------------*/

/* This preliminary study was looking at coincidence likelihoods to try
 * and sort events from non-events.  */
void likelihoodStudy(Run* mcs1, Run* mcs2, TH2D* lNPE, TH2D* lmaxT) {
	
	// initialize variables
	int cMode = mcs1->getCoincMode();
	double tTele = mcs1->getPeSumWindow();
	double lambda = 1/(2000.0*NANOSECOND); // time constant is 2us
	
	// make sure we're not just doing a 1 step (or aborted run)
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	if(dagSteps.size() < 2) {
		return;
	}
	
	// Foreground use dag steps 2 and 3
	bool background = false; // hardcoded background boolean
	double startTime;
	double endTime;
	
	if (background == false) {
		// find PMT hits for 2nd and 3rd dips (foreground)
		startTime  = *(dagSteps.begin() + 1);
		endTime = *(dagSteps.end() - 1);
		// force short runs
		if (startTime > 650 || endTime > 850) {
			return;
		}
	} else {
		// Find PMT hits for 750s during the long hold (background)
		startTime = dagSteps.front() - 850;
		endTime = dagSteps.front() - 100;
		if (startTime < 450) {
			return;
		}
	}
	// Make sure we don't have some weird timing glitch
	if (startTime >= endTime) {
		return;
	} else {
		printf("Found!\n"); 
	}
			
	// initialize and check coincidence vectors 
	// Determine which dagger pair we want to use
	Run* dagMCS;
	if ((mcs1->getCoincMode() == 1) || (mcs1->getCoincMode() == 2)) {
		dagMCS = mcs1;
	} else if ((mcs1->getCoincMode() == 3) || (mcs1->getCoincMode() == 4)) {
		dagMCS = mcs2;
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	// Load counts
	std::vector<coinc_t> coinc = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[startTime, endTime](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
	);
	std::vector<input_t> dagCts = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1, c2, startTime, endTime](input_t x)->bool{return ((x.ch == c1 || x.ch == c2) && x.realtime > startTime && x.realtime < (endTime+1.0));}
	);
	if(coinc.size() < 3) {
		return;
	}
	
	// Check to make sure beam doesn't turn on during backgrounds
	// While here check background rate, maybe?
	if((background == true)) {
		std::vector<input_t> oldCts = mcs1->getCounts(
			[](input_t x)->input_t{return x;},
			[startTime, endTime](input_t x)->bool{return x.ch == 3 && x.realtime > startTime && x.realtime < endTime;}
		);
		if ((double)oldCts.size() / (endTime - startTime) > 20.0) { // Arbitrarily say a 20hz rate is beam on!
			printf("Beam on during hold!\n");
			return;
		}	
	}
	
	// Load looping variables
	int phsA = 0;
	int phsB = 0;
	double tMax = 0.0;
	double like = 0.0;
	
	auto cIt = coinc.begin();
	for(auto dIt = dagCts.begin()+1; dIt < dagCts.end(); dIt++) {
		
		// Check that (a) the dagger is later than coincidence
		// (b) the dagger and previous dagger event are inside the coinc window
		// (c) we haven't hit the next coincidence event yet
		if ((dIt->time > cIt->time) && ((dIt->time - (dIt-1)->time) < tTele) && (dIt->time < (cIt+1)->time)) {
			
			// Add to Max time
			tMax = (double)(dIt->time- cIt->time)/0.8;
			// Increment dagger counter
			if((dIt->ch == 1) || (dIt->ch == 14)) {
				phsA += 1;
			} else {
				phsB += 1;
			}
			// Increment likelihood
			//like += log(lambda*exp(-lambda*(dIt->realtime - cIt->realtime)));
			like -= lambda*((double)(dIt->realtime - cIt->realtime));//*0.8*NANOSECOND;
			//like -= lambda*(dIt->realtime - cIt->realtime);
		
		// Dagger event before coincidence -- it's noise, just clear it and move on	
		} else if (dIt->time < cIt->time) {
			// clear data 
			phsA = 0;
			phsB = 0;
			tMax = 0.0;
			like = 0.0;
		// dagger event  no longer in coincidence! write out!
		} else if (((dIt->time - (dIt-1)->time) > tTele) && (dIt->time > cIt->time)){
		 		
			// make sure we can increment 
			if (cIt < coinc.end()) {
				cIt++;
			} else { break; }
			
			like += ((double)(phsA+phsB))*log(lambda);
			
			// Fill likelihood histograms
			lNPE->Fill((double)(phsA+phsB),like);
			lmaxT->Fill(tMax, like);
			if (like < -500.0) {
				printf("%f\n",like);
			}
			
			// clear data 
			phsA = 0;
			phsB = 0;
			tMax = 0.0;
			like = 0.0;	
		} /*else { // Weird things happening, fill and clear but don't move to next coinc.
			like += ((double)(phsA+phsB))*log(lambda);
			// Fill likelihood histograms
			lNPE->Fill((double)(phsA+phsB),like);
			lmaxT->Fill(tMax, like);
			if (like < -500.0) {
				printf("%f\n",like);
			}
			
			// clear data 
			phsA = 0;
			phsB = 0;
			tMax = 0.0;
			like = 0.0;	
		}*/
			
	}
}

/* 2D coincidence structure scan (coincidence length and photons). */
int multiCoincCheck(Run* mcs1, Run* mcs2, TH2D* nPhH, TH2D* coinL, TH2D* phByL) {
	
	// initialize variables
	int cMode = mcs1->getCoincMode();
	double fillEnd = getFillEnd(mcs1,mcs2);
	
	// make sure we're not just doing a 1 step (or aborted run)
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	if(dagSteps.size() < 2) {
		return 0;
	}
	
	double startTime = *(dagSteps.begin());
	double endTime   = *(dagSteps.end() - 1);

	double slop = 2.0; // hard coding slop of 2.0s
	double holdT = startTime - (fillEnd + 50.0); // 50s clean
	
	// hardcode in 20s run
	if (((holdT + slop) < 20.0) || ((holdT - slop) > 20.0)) {
		printf("Skipping Run -- not a 20s run!\n");
		return 0;
	}
	
	double lenUnld = (endTime - startTime);
	if (((lenUnld + slop) < 260.0) || ((lenUnld - slop) > 260.0)) {
		printf("Skipping run -- length of unload + bkgs not 260s!\n");
		return 0;
	}
			
	// initialize and check coincidence vectors
	Run* dagMCS;
	if ((mcs1->getCoincMode() % 2) == 1) {
		dagMCS = mcs1; // Low Threshold
	} else {
		dagMCS = mcs2; // High Threshold
	}
		
	// Load counts
	std::vector<coinc_t> coinc = dagMCS->getCoincCounts(
		[startTime](coinc_t x)->coinc_t{coinc_t y = x; y.realtime -= startTime; return y;},
		[startTime, endTime](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
	);
	
	if (coinc.size() <= 0) {
		printf("Stopping Run -- No coincidence counts found!\n");
		return 0;
	}
		
	// Fill Root Histogram
	for (auto cIt = coinc.begin(); cIt < coinc.end(); cIt++) {
		
		nPhH->Fill((10.0*cIt->realtime),((cIt->pmt1)+(cIt->pmt2)));
		coinL->Fill((10.0*cIt->realtime),(cIt->length/(10*NANOSECOND)));
		phByL->Fill((cIt->pmt1)+(cIt->pmt2),cIt->length/NANOSECOND);
	}
	
	printf("Loaded %lu coincidences!\n", coinc.size());
	return 1;
}


/* Object that saves the hits in one of our PMTs in a new photon tree */
//void savePMTHits(Run* mcs1, Run* mcs2, TH2D* arrivalTime, TTree* phsTree) {
void savePMTHits(Run* mcs1, Run* mcs2, TH1D* aTime1,TH1D* aTime2,TH1D* aTime3,TH1D* aTime4, TH2I* phsHist, int background) {
	// initialize variables and our new tree
	
	//phsTree->SetBranchAddress("a", (void *)&phsA);
	//phsTree->SetBranchAddress("b", (void *)&phsB);
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	// make sure we're not just doing a 1 step (or aborted run)
	if(dagSteps.size() < 2) {
		return;
	}
	
	// Foreground use dag steps 2 and 3
	//bool background = true;
	double startTime;
	double endTime;
	
	switch (background) { // This is the background summer thingy
		case 0: { // Background
			double bkgStart = *(dagSteps.end()-1) - 40;
			if (dagSteps.size() > 3) {
				bkgStart = mcs1->getTagBitEvt(1<<4, *(dagSteps.end()-2), 1) + 5.0;// 1<<4 is Cat Door, add extra 5 seconds to be safe
			}
			// and export our usual runs
			startTime = bkgStart;
			endTime   = *(dagSteps.end()-1);	
			// Only take long Run backgrounds
			if (startTime < 1000) {
				return;
			}
			printf("Bkg: %f, %f\n", startTime, endTime);
		} break;
		case 1: { // Production
			// find PMT hits for 2nd and 3rd dips (foreground)
			startTime  = *(dagSteps.begin() + 1);
			endTime = *(dagSteps.end() - 1);
			// force Short runs
			if (startTime < 650 || endTime < 850) {
				return;
			}
			printf("Prod.: %f, %f\n", startTime, endTime);
		} break;
		case 2: { // Hold
			// Find PMT hits for 250s during the long hold (background)
			startTime = dagSteps.front() - 350;
			endTime   = dagSteps.front() - 100;
			if (startTime < 450 || endTime <= startTime) {
				return;
			}
			
			double GV = getFillEnd(mcs1,mcs2);
			startTime = GV + 50; // New test -- first 250s of hold
			endTime   = GV + 300;
			/*startTime = dagSteps.front() - 20;
			endTime = dagSteps.front();
			if (startTime < 1000) {
				return;
			}*/
			printf("Hold: %f, %f\n",startTime,endTime);
		} break;
		case 3: { // End of Unload
			double countEnd = *(dagSteps.end()-1) - 50;
			if (dagSteps.size() > 3) {
				countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1); // TD move at end of run
			}
			if (dagSteps.front() < 1000) { // only short runs
				return;
			}
			startTime = countEnd - 50; // Go for last 50s
			endTime   = countEnd;
			printf("Unl: %f, %f\n",startTime,endTime);
		} break;
	}
	// Make sure we don't have some weird timing glitch
	if (endTime <= startTime) {
		printf("Start Time Glitch!");
		return;
	} else {
		printf("Found!\n"); 
	}
		
	// Load our dagger MCS depending on threshold
	int cMode = mcs1->getCoincMode();
	Run* dagMCS;
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	
	// Initialize and check coincidence vectors 
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	std::vector<coinc_t> coinc = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[startTime, endTime](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
	);
	std::vector<input_t> dagCts = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1,c2,startTime, endTime](input_t x)->bool{return ((x.ch == c1 || x.ch == c2) && x.realtime > startTime && x.realtime < (endTime+1.0));}
	);
	
	if(coinc.size() < 5) {
		return;
	}
	// Check to make sure beam doesn't turn on during backgrounds
	// While here check background rate, maybe?
	//if((background )) {
	std::vector<input_t> oldCts = mcs1->getCounts(
		[](input_t x)->input_t{return x;},
		[startTime, endTime](input_t x)->bool{return x.ch == 3 && x.realtime > startTime && x.realtime < endTime;}
	);
	if ((double)oldCts.size() / (endTime - startTime) > 20.0) { // Arbitrarily say a 20hz rate is beam on!
		printf("Beam on during hold!\n");
		return;
	}
	//}
	
	// parse the coincidence vector to get only coincidences 40us apart.
	// this is so we don't see photon tails, and thus avoid pileup.
	// WINDOW is a global of 40 us 
	std::vector<coinc_t> coincTimes;
	for(auto ii = coinc.begin() + 1; ii < coinc.end() - 1; ii++) {
		if((ii->realtime - (ii-1)->realtime) > WINDOW * NANOSECOND 
		   && ((ii+1)->realtime - ii->realtime) > WINDOW * NANOSECOND) {
			coincTimes.push_back(*ii);
		}
	}
	// doing a double-loopish thing here. 
	// choose the first coincidence we have
	auto cIt = coincTimes.begin();
	
	// Initialize pulse height structures
	int type;
	int phsA = 0;
	int phsB = 0;
	
	// Loop through all dagger counts
	int probA  = 0;
	//int echoPh = 0;
	int totalPh = 0;
	int totalC  = 0;
	for(auto dIt = dagCts.begin(); dIt < dagCts.end()-1; dIt++) {
		// If the dagger count is outside our window, it's not one of the coincidences we care about.
		// In this case reset! But first put our phs data into a previous histogram
		while((dIt->realtime - cIt->realtime) > WINDOW * NANOSECOND && cIt < coincTimes.end()) {
			cIt++;
			// Here we will fill the tree if there's a count
			if(phsA > 0 && phsB > 0) {
				phsHist->Fill(phsA,phsB);
				totalC += 1;
			}
			phsA = 0;
			phsB = 0;
			type = 0;
		}
		// Should have broken out of the dagger window loop. Now we fill photon height spectrum data
		if(dIt->realtime >= cIt->realtime && (dIt->realtime - cIt->realtime) <= WINDOW * NANOSECOND) {
			// Load the appropriate channel
			if((dIt->ch == 1) || (dIt->ch == 14)) {
				phsA += 1;
				totalPh+=1;
			}
			else {
				phsB += 1;
				totalPh+=1;
			}
			// This was a test for "echo" pulses. 
			//if (((dIt->realtime - cIt->realtime) <= 140.0 * NANOSECOND) &&  ((dIt->realtime - cIt->realtime) >= 120.0 * NANOSECOND)){
			//	echoPh +=1;
			//}
						
			// Figure out which PMT triggers first
			if(phsA == 1 && phsB == 0) { // PMTA triggers first
				type = 1;
				probA += 1;
			} 	else if(phsA == 0 && phsB == 1) {// PMTB triggers first
				type = 2;
			} 
			
			// Fill the event depending on what type of event it is.
			// ROOT hates everything, so I said fuck it and put in 4 different histograms.
			if(((dIt->ch == 1) || (dIt->ch == 14)) && type == 1 && phsA > 1) {aTime1->Fill(dIt->time - cIt->time);} 
			else if(((dIt->ch == 1) || (dIt->ch == 14)) && type == 2) {aTime2->Fill(dIt->time - cIt->time);} 
			else if(((dIt->ch == 2) || (dIt->ch == 15)) && type == 1) {aTime3->Fill(dIt->time - cIt->time);} 
			else if(((dIt->ch == 2) || (dIt->ch == 15)) && type == 2 && phsB > 1) {aTime4->Fill(dIt->time - cIt->time);}
		}
	}
	
	if (background) { // In background events, we'll write out PMT rate data
		
		FILE* outfile;
		outfile = fopen("PMTRates.csv", "a"); // Make sure we're cd'd into the right location
		
		std::vector<input_t> pmt1Cts = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[c1,startTime, endTime](input_t x)->bool{return (x.ch == c1 && x.realtime > startTime && x.realtime < endTime);}
		);
		std::vector<input_t> pmt2Cts = dagMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[c2,startTime, endTime](input_t x)->bool{return (x.ch == c2 && x.realtime > startTime && x.realtime < endTime);}
		);
		
		// These are the things I'm outputting to combine into a saved data file
		double pmt1Rate = pmt1Cts.size()/(endTime-startTime);
		double pmt2Rate = pmt2Cts.size()/(endTime-startTime);
		double coincRate = coincTimes.size()/(endTime-startTime);
		double pmtAProb  = ((double)probA) / ((double)totalC);
		
		//printf("Coinc Rate = %f\n", (double)coincTimes.size()/(endTime-startTime)); 
		fprintf(outfile,"%d,%f,%f,%f,%f\n",mcs1->getRunNo(),pmt1Rate,pmt2Rate,coincRate,pmtAProb);
		fclose(outfile);
	}
}

/* 2-dimensional pulse height vs. coincidence length of histogram producer.
 * Makes a "normal" and an "asymmetry" pulse height structure. */
void savePMTHitsTiming(Run* mcs1, Run* mcs2, TH2D* hist2D, TH2D* asym2D, int background) {
	// initialize variables and our new tree

	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	// make sure we're not just doing a 1 step (or aborted run)
	if(dagSteps.size() < 2) {
		return;
	}
	
	// Find the timing to use for our events based on "background" int
	bool SvsL = true; // Swap this to true for short
	double startTime;
	double endTime;
	switch (background) { // This is the background summer thingy
		case 0: { // Background
			double bkgStart = *(dagSteps.end()-1) - 40;
			if (dagSteps.size() > 3) {
				bkgStart = mcs1->getTagBitEvt(1<<4, *(dagSteps.end()-2), 1) + 5.0;// 1<<4 is Cat Door, add extra 5 seconds to be safe
			}
			// and export our usual runs
			startTime = bkgStart;
			endTime   = *(dagSteps.end()-1);
			if (SvsL) { // Short
				startTime = 10.; // Cheating -- Daytime Background
				endTime   = 250.;
				//if (startTime > 1000) { return;}
			} else { // Long
				if (startTime < 1000) {	return;}
			}
			//----------------------------------------------------------
			
			//----------------------------------------------------------
			printf("Bkg: %f, %f\n", startTime, endTime);
		} break;
		case 1: { // Production
			// find PMT hits for 2nd and 3rd dips (foreground)
			startTime  = *(dagSteps.begin() + 1);
			endTime = *(dagSteps.end() - 1);
			
			if (SvsL) { // Short
				startTime = 260.;
				endTime   = 500.;
				//if (startTime > 650 || endTime > 850) {	return;}
			} else { // Long
				if (startTime < 650 || endTime < 850) {	return;}
			}
			// force long runs
			
			printf("Prod.: %f, %f\n", startTime, endTime);
		} break;
		case 2: { // Hold
			// Find PMT hits for 250s during the long hold (background)
			if (SvsL) { // Short
				//startTime = dagSteps.front() - 20;
				//endTime = dagSteps.front();
				//if (startTime > 1000) {	return; }
				startTime = 510.;
				endTime   = 750.;
			} else { // Long
				startTime = dagSteps.front() - 350;
				endTime   = dagSteps.front() - 100;
				if (startTime < 450 || endTime <= startTime) { return; }
			}
			printf("Hold: %f, %f\n",startTime,endTime);
		} break;
		case 3: { // End of Unload
			double countEnd = *(dagSteps.end()-1) - 50;
			if (SvsL) { // Short
				/*if (dagSteps.size() > 3) {
					countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1); // TD move at end of run
				} 
				if (dagSteps.front() > 1000) { // only short runs
					return;
				}*/
				startTime = 760.;
				endTime   = 1000.;
			} else { // Long
				if (dagSteps.size() > 3) {
					countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1); // TD move at end of run
				} 
				if (dagSteps.front() < 1000) { // only long runs
					return;
				}
			}
			startTime = countEnd - 50; // Go for last 50s
			endTime   = countEnd;
			printf("Unl: %f, %f\n",startTime,endTime);
		} break;
		case 4: { // Background
			double bkgStart = *(dagSteps.end()-1) - 40;
			if (dagSteps.size() > 3) {
				bkgStart = mcs1->getTagBitEvt(1<<4, *(dagSteps.end()-2), 1) + 5.0;// 1<<4 is Cat Door, add extra 5 seconds to be safe
			}
			// and export our usual runs
			startTime = bkgStart;
			endTime   = *(dagSteps.end()-1);
			if (SvsL) { // Short
				if (startTime > 1000) { return;}
			} else { // Long
				if (startTime < 1000) {	return;}
			}
			//----------------------------------------------------------
			
			//----------------------------------------------------------
			printf("Bkg: %f, %f\n", startTime, endTime);
		} break;
		case 5: { // Production
			// find PMT hits for 2nd and 3rd dips (foreground)
			startTime  = *(dagSteps.begin() + 1);
			endTime = *(dagSteps.end() - 1);
			
			if (SvsL) { // Short
				if (startTime > 650 || endTime > 850) {	return;}
			} else { // Long
				if (startTime < 650 || endTime < 850) {	return;}
			}
			// force long runs
			
			printf("Prod.: %f, %f\n", startTime, endTime);
		} break;
		case 6: { // Hold
			// Find PMT hits for 250s during the long hold (background)
			if (SvsL) { // Short
				startTime = dagSteps.front() - 20;
				endTime = dagSteps.front();
				if (startTime > 1000) {	return; }
			} else { // Long
				startTime = dagSteps.front() - 350;
				endTime   = dagSteps.front() - 100;
				if (startTime < 450 || endTime <= startTime) { return; }
			}
			printf("Hold: %f, %f\n",startTime,endTime);
		} break;
		case 7: { // End of Unload
			double countEnd = *(dagSteps.end()-1) - 50;
			if (SvsL) { // Short
				if (dagSteps.size() > 3) {
					countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1); // TD move at end of run
				} 
				if (dagSteps.front() > 1000) { // only short runs
					return;
				}
			} else { // Long
				if (dagSteps.size() > 3) {
					countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1); // TD move at end of run
				} 
				if (dagSteps.front() < 1000) { // only long runs
					return;
				}
			}
			startTime = countEnd - 50; // Go for last 50s
			endTime   = countEnd;
			printf("Unl: %f, %f\n",startTime,endTime);
		} break;
		
	}
	// Make sure we don't have some weird timing glitch
	if (endTime <= startTime) {
		printf("Start Time Glitch!");
		return;
	} else {
		printf("Found!\n"); 
	}
		
	// Load our dagger MCS depending on threshold
	int cMode = mcs1->getCoincMode();
	Run* dagMCS;
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	
	// Going to do some SKETCH for a single PMT telescoping window
	
	// Initialize and check coincidence vectors 
	std::vector<coinc_t> coinc;
	
	if (cMode < 5) { // Fixed/Telescoping
		coinc = dagMCS->getCoincCounts(
			[](coinc_t x)->coinc_t{return x;},
			[startTime, endTime](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
		);
	}
	if (5 <= cMode < 9) { // Going to do 5/6 as PMT 1 and 7/8 as PMT 2
		// Integrated window
		double coinWindow = dagMCS->getCoincWindow(); // Initial window
		double PEWindow   = dagMCS->getPeSumWindow(); // Telescope window
		int PESum = ceil(dagMCS->getPeSum() / 2.0);  // Half the photon thr.
		if ((cMode==5) || (cMode==6)) { // PMT 1
			// Initialize and check coincidence vectors 
			int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
			std::vector<input_t> cts = dagMCS->getCounts(
				[](input_t x)->input_t{return x;},
				[c1,startTime, endTime](input_t x)->bool{return (x.ch == c1 && x.realtime > startTime && x.realtime < endTime);}
			);
			coinc = getSelfCoincs(cts, coinWindow, PEWindow, PESum);
		} else if ((cMode==7) || (cMode==8)) { // PMT 2
			// Initialize and check coincidence vectors 
			int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
			std::vector<input_t> cts = dagMCS->getCounts(
				[](input_t x)->input_t{return x;},
				[c2,startTime, endTime](input_t x)->bool{return (x.ch == c2 && x.realtime > startTime && x.realtime < endTime);}
			);
			coinc = getSelfCoincs(cts, coinWindow, PEWindow, PESum);
		}
	}
	
	
	if(coinc.size() < 5) {
		return;
	}
	// Check to make sure beam doesn't turn on during backgrounds
	// While here check background rate, maybe?
	std::vector<input_t> oldCts = mcs1->getCounts(
		[](input_t x)->input_t{return x;},
		[startTime, endTime](input_t x)->bool{return x.ch == 3 && x.realtime > startTime && x.realtime < endTime;}
	);
	if ((double)oldCts.size() / (endTime - startTime) > 20.0) { // Arbitrarily say a 20hz rate is beam on!
		printf("Beam on during hold!\n");
		return;
	}
		
	// parse the coincidence vector to get only coincidences 40us apart.
	// this is so we don't see photon tails, and thus avoid pileup.
	// WINDOW is a global of 40 us 
	std::vector<coinc_t> coincTimes;
	for(auto ii = coinc.begin() + 1; ii < coinc.end() - 1; ii++) {
		if((ii->realtime - (ii-1)->realtime) > WINDOW * NANOSECOND 
		   && ((ii+1)->realtime - ii->realtime) > WINDOW * NANOSECOND) {
			coincTimes.push_back(*ii);
		}
	}
	
	// Looping through coincidences
	for (auto cIt = coincTimes.begin(); cIt < coincTimes.end()-1; cIt++) {
		
		int total  = cIt->pmt1 + cIt->pmt2; 
		double len = cIt->length / NANOSECOND;
		//double ratio = (double)cIt->pmt1 / ((double)total);
		int fast = cIt->prompt;
		//printf("%d,%d,%f\n",total,fast,len);
		// Fill the two histograms
		hist2D->Fill(total,len);
		asym2D->Fill(total,fast);		
	}
}

void savePMTHits2D(Run* mcs1, Run* mcs2, TH2D* aTime1,TH2D* aTime2,TH2D* aTime3,TH2D* aTime4, bool background) {
	// initialize variables and our new tree
	int type;
	int cMode = mcs1->getCoincMode();
	int phsA = 0;
	int phsB = 0;
	
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	// make sure we're not just doing a 1 step (or aborted run)
	if(dagSteps.size() < 2) {
		return;
	}
	
	// Foreground use dag steps 2 and 3
	//bool background = true;
	double startTime;
	double endTime;
	
	if (background == false) {
		// find PMT hits for 2nd and 3rd dips (foreground)
		startTime  = *(dagSteps.begin() + 1);
		endTime = *(dagSteps.end() - 1);
		// force short runs
		if (startTime > 650 || endTime > 850) {
			return;
		}
		printf("%f, %f\n", startTime, endTime);
	} else {
	
	// Find PMT hits for 750s during the long hold (background)
		startTime = dagSteps.front() - 850;
		endTime = dagSteps.front() - 100;
		if (startTime < 450) {
			return;
		}
	}
	// Make sure we don't have some weird timing glitch
	if (startTime >= endTime) {
		return;
	} else {
		printf("Found!\n"); 
	}
		
	// initialize and check coincidence vectors
	Run* dagMCS;
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
		
	// Load counts
	std::vector<coinc_t> coinc = dagMCS->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[startTime, endTime](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
	);
	
	std::vector<input_t> dagCts = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1,c2, startTime, endTime](input_t x)->bool{return ((x.ch == c1 || x.ch == c2) && x.realtime > startTime && x.realtime < (endTime+1.0));}
	);
	
	if(coinc.size() < 3 || dagCts.size() < 10) {
		printf("Stopping Run -- No coincidence counts found!\n");
		printf("%lu,%lu\n",coinc.size(),dagCts.size());
		return;
	}
	printf("%lu\n",coinc.size());
	// Check to make sure beam doesn't turn on during backgrounds
	// While here check background rate, maybe?
	if((background == true)) {
		std::vector<input_t> oldCts = mcs1->getCounts(
			[](input_t x)->input_t{return x;},
			[startTime, endTime](input_t x)->bool{return x.ch == 3 && x.realtime > startTime && x.realtime < endTime;}
		);
		if ((double)oldCts.size() / (endTime - startTime) > 20.0) { // Arbitrarily say a 20hz rate is beam on!
			printf("Beam on during hold!\n");
			return;
		}
	}
	
	// parse the coincidence vector to get only coincidences 40us apart.
	// this is so we don't see photon tails, and thus avoid pileup.
	//double WINDOW = 40000; // 40us
	std::vector<coinc_t> coincTimes;
	for(auto ii = coinc.begin() + 1; ii < coinc.end() - 1; ii++) {
		if((ii->realtime - (ii-1)->realtime) > WINDOW * NANOSECOND 
		   && ((ii+1)->realtime - ii->realtime) > WINDOW * NANOSECOND) {
			coincTimes.push_back(*ii);
		}
	}
	
	// doing a double-loopish thing here. 
	// choose the first coincidence we have
	auto cIt = coincTimes.begin();
	// Loop through all dagger counts
	for(auto dIt = dagCts.begin(); dIt < dagCts.end()-1; dIt++) {
		// If the dagger count is outside our window, it's not one of the coincidences we care about.
		// In this case reset! But first put our phs data into a previous histogram
		while(((dIt->realtime - cIt->realtime) > WINDOW * NANOSECOND) && (cIt < coincTimes.end())) {
			cIt++;
			// Here we will fill the tree if there's a count
			if(phsA > 0 && phsB > 0) {
				//phsHist->Fill(phsA,phsB);
				
				//phsTree->Fill();
			}
			phsA = 0;
			phsB = 0;
			type = 0;
		}
		// Should have broken out of the dagger window loop. Now we fill photon height spectrum data
		if ((dIt->realtime >= cIt->realtime) && ((dIt->realtime - cIt->realtime) <= WINDOW * NANOSECOND)) {
			// Load the appropriate channel
			if(dIt->ch == c1) {
				phsA += 1;
			}
			else if (dIt->ch == c2) {
				phsB += 1;
			}
			// Figure out which PMT triggers first, that's what "type" means 
			if(phsA == 1 && phsB == 0) {type = 1;} // PMTA triggers first
			else if(phsA == 0 && phsB == 1) {type = 2;} // PMTB triggers first
			
			// Fill the event depending on what type of event it is.
			// ROOT hates everything, so I said fuck it and put in 4 different histograms
			//if ((cIt->avg1 < 17*NANOSECOND) || (cIt->avg2 < 17*NANOSECOND)) {
			if ((cIt->avg1 > 0 && cIt->avg1 < 20 * NANOSECOND) || (cIt->avg2 > 0 && cIt->avg2 < 20*NANOSECOND)) {
				continue; } else {
			printf("%d %f %d %f %f\n",cIt->pmt1,cIt->avg1 / NANOSECOND, cIt->pmt2,cIt->avg2 / NANOSECOND, cIt->length / NANOSECOND);
			
			if(dIt->ch == c1 && (type == 1 && phsA > 1)) {aTime1->Fill(dIt->time - cIt->time, phsA);} 
			else if(dIt->ch == c1 && type == 2) {aTime2->Fill(dIt->time - cIt->time, phsA);} 
			else if(dIt->ch == c2 && type == 1) {aTime3->Fill(dIt->time - cIt->time, phsB);} 
			else if(dIt->ch == c2 && (type == 2 && phsB > 1)) {aTime4->Fill(dIt->time - cIt->time, phsB);}
			
		}
		}
	}
}

/* Create Histogram of Next PMT hit */
void singlePMTHits(Run* mcs1, Run* mcs2, TH1D* dT1, TH1D* dT2) {

	// make sure we're not just doing a 1 step (or aborted run)
	int cMode = mcs1->getCoincMode();
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	if(dagSteps.size() < 2) {
		return;
	}
		
	// Find PMT hits for 750s during the long hold (background)
	//double startTime = dagSteps.front() - 850;
	//double endTime = dagSteps.front() - 100;
	double startTime  = *(dagSteps.begin() + 1);
	double endTime = *(dagSteps.end() - 1);
	if (startTime < 450) {
		return;
	}
	
	// Make sure we don't have some weird timing glitch
	if (startTime >= endTime) {
		return;
	} else {
		printf("Found!\n"); 
	}
	
	// Take all photon counts on PMT1 + PMT2
	// initialize and check coincidence vectors
	Run* dagMCS;
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	std::vector<input_t> d1Cts = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1, startTime, endTime](input_t x)->bool{return ((x.ch == c1) && x.realtime > startTime && x.realtime < endTime);}
	);
	std::vector<input_t> d2Cts = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2, startTime, endTime](input_t x)->bool{return ((x.ch == c2) && x.realtime > startTime && x.realtime < endTime);}
	);
	
	// Loop through with detector iterators and determine our separations
	double sep = 0.0;
	for (auto d1It = d1Cts.begin(); d1It < d1Cts.end(); d1It ++) {
		
		sep = ((d1It+1)->time) - (d1It->time);
		if (sep > 0.0) {
			dT1->Fill(sep);
		}
	}

	for (auto d2It = d2Cts.begin(); d2It < d2Cts.end(); d2It ++) {	
		sep = ((d2It+1)->time) - (d2It->time);
		if (sep > 0.0) {
			dT2->Fill(sep);
		}
	}	
	
}

/* Look at PMT structure during the fill to make sure there's not 
 * double pulsing or anything going on there.
 * 
 * There's some sketchy "coincidence mode" encoding as the channel.  */

void fillDetHitsMoving(Run* mcs1, Run* mcs2, TH1D* dT) {
//void fillDetHitsMoving(Run* mcs1, Run* mcs2, TH1D* dT,TH2D* hist2D) {
	
	// I'm encoding channel in the coincidence mode.
	int chan = mcs1->getCoincMode();
	
	// Find Timing from first H-GX pulse to end of fill
	double startTime = mcs2->getTagBitEvt(1<<1, 0.0, 1);
	double endTime = getFillEnd(mcs1, mcs2);
	
	if (((startTime < 0.0) || (endTime < 0.0)) || (startTime >= endTime)) {
		printf("Error: Cannot find hits or fill end!\n");
		return;
	}
	
	double multiplier;
	if ((chan == 1 || chan == 2) || (chan == 6) || (chan == 7 || chan == 8)) {
		multiplier = 1.0;
	} else {
		multiplier = 800.0;
	}
		
	Run* detMCS;
	int detChan;
	if (chan <= 5) {
		detMCS = mcs1;	
		detChan = chan + detMCS->getMCSOff();
	} else if ((5 < chan) && (chan <= 10)) {
		detMCS = mcs2;
		detChan = (chan-5) + detMCS->getMCSOff();
	}
	
	// Take all photon hits on our detector
	std::vector<input_t> detCts;
	detCts = detMCS->getCounts(
			[](input_t x)->input_t{return x;},
			[detChan, startTime, endTime](input_t x)->bool{return ((x.ch == detChan) && 
															(x.realtime > startTime && x.realtime < endTime));}
			);
	//std::sort(detCts.begin(),detCts.end(),[](input_t x, input_t y)->bool{return (x.realtime < y.realtime);});
	double sep = 0.0;
	double maxSep = 5000.0; // Max separation of 5000 MCS ticks (4 us)
	// Loop through photon counts
	for (auto cIt = detCts.begin(); cIt < detCts.end(); cIt++) {
		for (auto dIt = cIt+1; dIt < detCts.end(); dIt++) {
			// Check the separation is between zero and our fixed maximum
			sep = ((dIt->time) - (cIt->time)) * multiplier; // Separation in MCS ticks
			if (sep > maxSep) { break; }
			
			if ((sep > 0.0) && (sep < maxSep)) {
				dT->Fill(sep);
			}
		}
	}
	
	//std::vector<coinc_t> coincTest = getSelfCoincs(detCts, 1000, 10000, 4);
	//printf("coincTest size: %lu\n", coincTest.size());
}

/* Create ROOT Histogram of all PMT hits within 1us after an initial hit */
void singlePMTHitsMoving(Run* mcs1, Run* mcs2, TH1D* dT1, TH1D* dT2) {
	
	// make sure we're not just doing a 1 step (or aborted run)
	int cMode = mcs1->getCoincMode();
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	if(dagSteps.size() < 2) {
		return;
	}
		
	// Find PMT hits for 750s during the long hold (background)
	double startTime = dagSteps.front() - 850;
	double endTime = dagSteps.front() - 100;
	//double startTime  = *(dagSteps.begin() + 1);
	//double endTime = *(dagSteps.end() - 1);
	if (startTime < 450) {
		return;
	}
	
	// Make sure we don't have some weird timing glitch
	if (startTime >= endTime) {
		return;
	} else {
		printf("Found!\n"); 
	}
	
	// Take all photon counts on PMT1 + PMT2
	// initialize and check coincidence vectors
	Run* dagMCS;
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	std::vector<input_t> d1Cts = mcs1->getCounts(
		[](input_t x)->input_t{return x;},
		[c1, startTime, endTime](input_t x)->bool{return ((x.ch == c1) && x.realtime > startTime && x.realtime < endTime);}
	);
	std::vector<input_t> d2Cts = mcs1->getCounts(
		[](input_t x)->input_t{return x;},
		[c2, startTime, endTime](input_t x)->bool{return ((x.ch == 2) && x.realtime > startTime && x.realtime < endTime);}
	);
		
	double sep = 0.0;
	double maxSep = 1000.0; // Max separation of 1000 nanoseconds
	// Loop through both PMT counters
	if (d1Cts.size() > 0) { // ignore if there's no counts
		for (auto cIt = d1Cts.begin(); cIt < d1Cts.end(); cIt++) {
			for (auto dIt = cIt; dIt <= d1Cts.end(); dIt++) {
				// Check the separation is between zero and our fixed maximum
				sep = ((dIt->time) - (cIt->time));
				if (sep > maxSep) { break; }
				
				if ((sep > 0.0) && (sep < maxSep)) {
					dT1->Fill(sep);
				}
				//if ((sep > 120.0 * NANOSECOND) && (sep < 160.0 * NANOSECOND)) {
					//dT2->Fill(cIt->realtime);
				//}		 
			}
		}
	}
	if (d2Cts.size() > 0) { // ignore if there's no counts
		for (auto cIt = d2Cts.begin(); cIt < d2Cts.end(); cIt++) {
			for (auto dIt = cIt; dIt <= d2Cts.end(); dIt++) {
				// Check the separation is between zero and our fixed maximum
				sep = ((dIt->time) - (cIt->time));
				if (sep > (maxSep)) { break; } 
				
				if ((sep > 0.0 && sep < maxSep)) {
					dT2->Fill(sep);
				}
				//if ((sep > 120.0 * NANOSECOND) && (sep < 160.0 * NANOSECOND)) {
					//dT2->Fill(cIt->realtime);
					//printf("%e\n",cIt->realtime);
				//}
				
			}
		}
	}
}

/* This does nothing. */
void writePMTBalance(Run* mcs1, Run* mcs2) {
	
	// initialize variables and our new tree
	double pmtAProb;
	int type;
	int phsA = 0;
	int phsB = 0;
	
	// Load info from Run*
	int cMode = mcs1->getCoincMode();
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	// make sure we're not just doing a 1 step (or aborted run)
	if(dagSteps.size() < 2) {
		return;
	}
	
	// For this run we only care about the foreground events -- want
	double startTime;
	double endTime;
	
	
	
	//#define BKGRATEA 560 // Singles background rates PMT1
	//#define BKGRATEB 545 // Singles background rates PMT2
	//#define COINCBKGRATE 4.89 // Background rates from coincidence
	//#define PMTAPROB 0.584975 // Need to calculate this still!

	
	
	
	
	// Foreground use dag steps 2 and 3
	bool background = false;
	/*
	if (background == false) {
		// find PMT hits for 2nd and 3rd dips (foreground)
		startTime  = *(dagSteps.begin() + 1);
		endTime = *(dagSteps.end() - 1);
		// force short runs
		
		if (startTime > 650 || endTime > 850) {
			return;
		}
		printf("%f, %f\n", startTime, endTime);
	} else {
	
	// Find PMT hits for 750s during the long hold (background)
		startTime = dagSteps.front() - 850;
		endTime = dagSteps.front() - 100;
		if (startTime < 450) {
			return;
		}
	}
	// Make sure we don't have some weird timing glitch
	if (startTime >= endTime) {
		return;
	} else {
		printf("Found!\n"); 
	}
	
	// initialize and check coincidence vectors 
	std::vector<input_t> coinc;
	std::vector<input_t> dagCts;
	if ((cMode == 1) || (cMode == 2)) {
		coinc = mcs1->getCoincCounts(
			[](input_t x)->input_t{return x;},
			[startTime, endTime](input_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
		);
		dagCts = mcs1->getCounts(
			[](input_t x)->input_t{return x;},
			[startTime, endTime](input_t x)->bool{return ((x.ch == 1 || x.ch == 2) && x.realtime > startTime && x.realtime < (endTime+1.0));}
		);
	} else if ((cMode == 3) || (cMode == 4)) { // High threshold
		coinc = mcs2->getCoincCounts(
			[](input_t x)->input_t{return x;},
			[startTime, endTime](input_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
		);
		dagCts = mcs2->getCounts(
			[](input_t x)->input_t{return x;},
			[startTime, endTime](input_t x)->bool{return ((x.ch == 14 || x.ch == 15) && x.realtime > startTime && x.realtime < (endTime+1.0));}
		);
	}
	
	if(coinc.size() < 20) { // For some reason there's a 16 count PMT background (13249) where that PMT should be off
		return;
	}
	// Check to make sure beam doesn't turn on during backgrounds
	// While here check background rate, maybe?
	if((background == true)) {
		std::vector<input_t> oldCts = mcs1->getCounts(
			[](input_t x)->input_t{return x;},
			[startTime, endTime](input_t x)->bool{return x.ch == 3 && x.realtime > startTime && x.realtime < endTime;}
		);
		if ((double)oldCts.size() / (endTime - startTime) > 20.0) { // Arbitrarily say a 20hz rate is beam on!
			printf("Beam on during hold!\n");
			return;
		}
		
	}
	
	// parse the coincidence vector to get only coincidences 40us apart.
	// this is so we don't see photon tails, and thus avoid pileup.
	//double WINDOW = 40000; // 40us
	std::vector<input_t> coincTimes;
	for(auto ii = coinc.begin() + 1; ii < coinc.end() - 1; ii++) {
		if((ii->realtime - (ii-1)->realtime) > WINDOW * NANOSECOND 
		   && ((ii+1)->realtime - ii->realtime) > WINDOW * NANOSECOND) {
			coincTimes.push_back(*ii);
		}
	}
	// doing a double-loopish thing here. 
	// choose the first coincidence we have
	auto cIt = coincTimes.begin();
	// Loop through all dagger counts
	int echoPh = 0;
	int totalPh = 0;
	for(auto dIt = dagCts.begin(); dIt < dagCts.end()-1; dIt++) {
		// If the dagger count is outside our window, it's not one of the coincidences we care about.
		// In this case reset! But first put our phs data into a previous histogram
		while((dIt->realtime - cIt->realtime) > WINDOW * NANOSECOND && cIt < coincTimes.end()) {
			cIt++;
			// Here we will fill the tree if there's a count
			if(phsA > 0 && phsB > 0) {
				phsHist->Fill(phsA,phsB);
				//phsTree->Fill();
			}
			phsA = 0;
			phsB = 0;
			type = 0;
		}
		// Should have broken out of the dagger window loop. Now we fill photon height spectrum data
		if(dIt->realtime >= cIt->realtime && (dIt->realtime - cIt->realtime) <= WINDOW * NANOSECOND) {
			// Load the appropriate channel
			if((dIt->ch == 1) || (dIt->ch == 14)) {
				phsA += 1;
				totalPh+=1;
			}
			else {
				phsB += 1;
				totalPh+=1;
			}
			if (((dIt->realtime - cIt->realtime) <= 140.0 * NANOSECOND) &&  ((dIt->realtime - cIt->realtime) >= 120.0 * NANOSECOND)){
				echoPh +=1;
			}
			
			
			// Figure out which PMT triggers first, that's what "type" means 
			if(phsA == 1 && phsB == 0) {type = 1;} // PMTA triggers first
			else if(phsA == 0 && phsB == 1) {type = 2;} // PMTB triggers first
			
			// Fill the event depending on what type of event it is.
			// ROOT hates everything, so I said fuck it and put in 4 different histograms.
			if(((dIt->ch == 1) || (dIt->ch == 14)) && type == 1 && phsA > 1) {aTime1->Fill(dIt->time - cIt->time);} 
			else if(((dIt->ch == 1) || (dIt->ch == 14)) && type == 2) {aTime2->Fill(dIt->time - cIt->time);} 
			else if(((dIt->ch == 2) || (dIt->ch == 15)) && type == 1) {aTime3->Fill(dIt->time - cIt->time);} 
			else if(((dIt->ch == 2) || (dIt->ch == 15)) && type == 2 && phsB > 1) {aTime4->Fill(dIt->time - cIt->time);}
		}
	}
*/	
}


/* Another background object -- this one for TTNE. 
 * Deprecated.*/
void bkgTTNE(Run* run, TH1D* ttneHist) {
	std::vector<input_t> cts = run->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return (x.ch == 1 || x.ch == 2) && x.realtime > 200.0 && x.realtime < 1200.0;}
	);
	std::vector<coinc_t> coinc = run->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[](coinc_t x)->bool{return x.realtime > 200.0 && x.realtime < 1200.0;}
	);
	auto cIt = coinc.begin();
	for(auto it = cts.begin(); it < cts.end()-1; it++) {
		// If we've passed the coincidence time look for the next one.
		while(cIt < coinc.end()-1 && floor(cIt->realtime) < floor(it->realtime)) { 
			cIt++;
		}
		if(floor(it->realtime) == floor(cIt->realtime)) {
			continue;
		}
		if(it->ch == 2) {
			continue;
		}
		auto nIt = it+1;
		while(nIt < cts.end() && nIt->ch != 2) {
			nIt++;
		}
		if(nIt < cts.end() && nIt->ch == 2 && floor(nIt->realtime) != floor(cIt->realtime)) {
			ttneHist->Fill(nIt->time - it->time);
		}
	}
}

/* Get ratios and rates on Ttne. Deprecated. */
void getTtne(Run* run, TH1D* ttneHist, int ch, double ts, double te) {
	std::vector<input_t> cts = run->getCounts(
		[](input_t x)->input_t{return x;},
		[ch, ts, te](input_t x)->bool{return x.ch == ch && x.realtime > ts && x.realtime < te;}
	);

	uint64_t sumDoub = 0;
	uint64_t sumReg = 0;
	for(auto it = cts.begin()+1; it < cts.end(); it++) {
		double diff = (it->realtime - (it-1)->realtime) / NANOSECOND;
		if(diff < 10000) {
			sumDoub += 1;
		}
		else {
			sumReg += 1;
		}
	}
	printf("Ratio - %d %f %lu %lu\n", run->getRunNo(), ((double)sumDoub)/((double)sumReg), sumDoub, sumReg);
	printf("Rate - %f\n", cts.size()/(te-ts));
}
