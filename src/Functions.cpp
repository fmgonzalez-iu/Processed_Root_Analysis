#include "../inc/Functions.hpp"
#include "../inc/ExpFill.hpp"
#include "TCanvas.h"
#include "TGraph.h"

/* define constants we need for later */
#define NANOSECOND .000000001
//#define WINDOW 40000
#define bkgMov50ns8pe 0.10666
#define bkgMov50ns1000ns6pe 0.3
#define synthbkg_50_500_2 0.0
#define TAUN 877.7

/*----------------------------------------------------------------------------
	Author: Nathan B. Callahan (?)
	Editor: Frank M. Gonzalez
	
	This file contains a bunch of functions that we use in other parts of
	the analysis scripts. 
------------------------------------------------------------------------------*/

/* Dead-time exists in the detector due to finite discriminator width.
 *  n = m/(1-m\tau) 
 * where n is the corrected rate, m is the recorded rate, and \tau is the fixed deadtime*/
 /*
double getDeadTimeCoinc(std::vector<coinc_t> &cts, double start, double end) {
	
	//int coincType = run->getCoincMode();

	double deadTimeCounts = 0.0;
	//std::vector<coinc_t> cts = run->getCoincCounts(
	//		[](coinc_t x)->coinc_t{return x;},
	//		[start, end](coinc_t x)->bool{return (x.realtime > start && x.realtime < end);});
			
	if(cts.size() < 1) {
		return deadTimeCounts;
	}

	int counts = 0;
	double totCoincL = 0.0;
	double dtS = 0.1; // Divide into 0.1s bins
	
	// Deadtime is assumed to be the average of the coincidence window length
	int time=floor(cts.front().realtime / dtS);
	//double peSumWindow = run->getPeSumWindow();
	for(auto ctsIterator = cts.begin(); ctsIterator < cts.end(); ctsIterator++) {
		if(floor((*ctsIterator).realtime / dtS) > time) {
			// n = m / (1 - m tau)
			// n = (cts / dt) / (1 - (cts / dt) (totTau / cts) )
			// n = (cts / dt) / (1 - totTau / dt) 
			// dtC = n * dt 
			// dtC = cts / ( 1 - totTau / dt)
			deadTimeCounts += (double)counts / (1.0 - totCoincL / dtS)  - (double)counts;
			counts = 0;
			totCoincL = 0.0;
			time=floor((*ctsIterator).realtime / dtS);
		}
		totCoincL += (*ctsIterator).length;
		counts++;
	}
	
	return deadTimeCounts;
}

double getDeadTimeSing(std::vector<input_t> &cts, double start, double end) {
	// Singles deadtime. Deadtime is a part of the input_t class
	// hardcoded in 16ns unless otherwise changed

	// deadtime counts defined here
	double deadTimeCounts = 0.0;
	if(cts.size() < 1) {
		return deadTimeCounts;
	}

	int counts = 0;
	double totDT = 0.0;
	double dtS = 0.1; // Divide into 0.1s bins
	
	// Deadtime is assumed to be the average of the coincidence window length
	int time=floor(cts.front().realtime / dtS);
	//double peSumWindow = run->getPeSumWindow();
	bool test = false;
	for(auto cIt = cts.begin(); cIt < cts.end(); cIt++) {
		
		//if ((*cIt).deadtime != 16*NANOSECOND) { test=true; printf("Deadtime: %f\nRate:%f\n",(*cIt).deadtime/NANOSECOND,(*cIt).rate);}
		if(floor((*cIt).realtime / dtS) > time) {
			//if (test) { printf("counts %d, totDT %f\n",counts,totDT); } 
			deadTimeCounts += (double)counts / (1.0 - totDT / dtS)  - (double)counts;
			counts = 0;
			totDT = 0.0;
			time = floor((*cIt).realtime / dtS);
			test = false;
		}
		totDT += (*cIt).deadtime; // Add the deadtime of each event
		counts++;
	}
	
	return deadTimeCounts;		
			
}

/* For singles, we have a deadtime due to hardware */
/*double getDeadTimeSing(std::vector<input_t> &cts) { 
	
	// Initialize, check that we have counts
	double deadTimeCts = 0.;
	if (cts.size() < 1) {
		return deadTimeCts;
	}
	
	double stepMean1 = cts1.size() > 0 ? std::accumulate(cts1.begin(), cts1.end(), 0.0, 
							[](double m, input_t x)->double{return m + x.realtime;}) / (double)cts1.size() : 0.0;
	
	for(auto ctsIterator = cts.begin(); ctsIterator < cts.end(); ctsIterator++) {
		if(floor((*ctsIterator).realtime / dtS) > time) {
			// n = m / (1 - m tau)
			// n = (cts / dt) / (1 - (cts / dt) (totTau / cts) )
			// n = (cts / dt) / (1 - totTau / dt) 
			// dtC = n * dt 
			// dtC = cts / ( 1 - totTau / dt)
			deadTimeCounts += (double)counts / (1.0 - totCoincL / dtS)  - (double)counts;
			counts = 0;
			totCoincL = 0.0;
			time=floor((*ctsIterator).realtime / dtS);
		}
		totCoincL += (*ctsIterator).length;
		counts++;
	}
	
	
	
	return deadTimeCts;
	
}*/

/* This is a single-exponential pileup */
/*double getPileUpTau(Run* run, std::vector<coinc_t> coin_vec, double bkgRate, double start, double end) {
	// Assume that all photons above "background" are from "coincidences"
	// that are uncounted by our algorithm. Also assume the form of coinc.
	// events are the sum of exponentials. Then, find:
	// \frac{x - y}{z} = \frac{N_m}{\Sum_{j}{\kappa_j \tau_j \left(1 - e^{-\delta t_m / \tau_j}\right)}}
	//
	// Here we're solving for the "lifetime" for scintillator glow, assuming single exponential

	double tau = 0.;
	if (coin_vec.size() == 0) { return tau; } // Need coincidences!

	// Load singles from data.
	// Have to treat PMT1 and PMT2 separately due to deadtime 
	int c1 = run->getCoinC1() + run->getMCSOff();
	std::vector<input_t> pmt1_vec = run->getCounts(
		[](input_t x)->input_t{return x;},
		[c1,start,end](coinc_t x)->bool{return (x.ch  == c1) 
										 && (x.realtime > start && x.realtime < end);}
	);
	double pmt1 = (double)pmt1_vec.size() + getDeadTimeSing(pmt1_vec, start, end);
	
	int c2 = run->getCoinC2() + run->getMCSOff();
	std::vector<input_t> pmt2_vec = run->getCounts(
		[](input_t x)->input_t{return x;},
		[c2,start,end](coinc_t x)->bool{return (x.ch  == c2) 
										 && (x.realtime > start && x.realtime < end);}
	);
	double pmt2 = (double)pmt2_vec.size() + getDeadTimeSing(pmt2_vec, start, end);
	
	// Previously defined coincidence data 
	// (so that this is generizable to pseudo-coincidences)
	double coin = (double)coin_vec.size() + getDeadTimeCoinc(coin_vec, start, end);
	
	// Also need to know average length and photons of coincidences:
	double dt_m = std::accumulate(coin_vec.begin(), coin_vec.end(), 0.0, [](double m, coinc_t x)->double{return m + (x.length);}) / (double)coin_vec.size();// Average length of coincidence
	double N_m  = std::accumulate(coin_vec.begin(), coin_vec.end(), 0.0, [](double m, coinc_t x)->double{return m + (double)(x.pmt1 + x.pmt2);}) / (double)coin_vec.size(); // Average number of photons in coincidences
		
	// Over the whole unload, want to know:
	double ph_f = pmt1 + pmt2; // Number of photons in foregrounds
	double ph_b = bkgRate * (end - start); // Number of photons in backgrounds
	
	double amp = (ph_f - ph_b) / coin; // Amplitude of average coinc. event
		
	// Exponential lifetime of coincidence glow
	tau = - dt_m /  log( 1. - (N_m / amp));
	
	return tau;
}

/* How many coincidences do we lose, given single tau */
/*double getPileUpCorr(std::vector<coinc_t> coinc, double tau, double amp, int PE) { 
	// With a given (single) tau, we need to deweight coincidences
	// based on their interarrival times. 
	//
	// Rather than doing a factorial
	
	
	double poisson = mu_i**pe * exp(mu_i) / pe!
	
}*/

double getACUp(Run* mcs1) {
	
	// 2 second timing slop (arbitrary)
	double slop = 2.0;
	// These events should happen once per run
	double acUp    = mcs1->getTagBitEvt(1<<0, slop, 1);
	
	return acUp;
}

/* Determine the end of fill, based on tagbits. Differs depending on the run. */
double getFillEnd(Run* mcs1, Run* mcs2) {

	double moveErr = 0.2; // slop in the timings
	double fillEnd;
	// hard coded in backup fillEnds
	double fillEnd_RH = 300.0;
	double fillEnd_noRH = 150.0; 
	
	// get the run number (useful for hardcoding in events -- e.g. tagbits off or something...)
	int runNo = mcs1->getRunNo();
	
	// Get some tagbit times. These should all be the same.
	// Priority: GV->TD->Dagger->HardCode
	double gvTime = mcs1->getTagBitEvt(1<<2, 10.0, 1);
	double tdTime = mcs1->getTagBitEvt(1<<4, 10.0, 1); // technically Cat Door. This is the first thing to move.
	double duTime = mcs2->getTagBitEvt(1<<0, 10.0, 1);
	
	// Easiest thing: GV and TD agree. This only works for production.
	if ((gvTime > 0) && (fabs(gvTime - tdTime) < moveErr)) {
		fillEnd = gvTime;
	// Next check: if TD and Dagger movement agree (for e.g. backgrounds)
	} else if ((tdTime > 0) && (fabs(tdTime - duTime) < moveErr)) {
		fillEnd = tdTime;
	// Now we go to hardcoded-ish tests. If there's any movement around a hard-coded set. 
	} else if ((fabs(gvTime - fillEnd_noRH) < moveErr) || (fabs(gvTime - fillEnd_RH) < moveErr)) {
		fillEnd = gvTime;
	} else if ((fabs(tdTime - fillEnd_noRH) < moveErr) || (fabs(tdTime - fillEnd_RH) < moveErr)) {
		fillEnd = tdTime;
	} else if ((fabs(duTime - fillEnd_noRH) < moveErr) || (fabs(duTime - fillEnd_RH) < moveErr)) {
		fillEnd = duTime;
	// Now we're just going to pure hardcoding based on run number. 
	// These run numbers are probably wrong for some of these.
	} else if (runNo < 9020 || ((runNo < 14711) && (runNo > 14722))) {
		printf("Using constant fillEnd (%fs) for this run !!!\n", fillEnd_noRH);
		fillEnd = fillEnd_noRH;
	} else {
		printf("Using constant fillEnd (%fs) for this run !!!\n", fillEnd_RH);
		fillEnd = fillEnd_RH;
	}
	
	return fillEnd;
}

/* Find the unload dagger steps  */
std::vector<double> dagDips(Run* mcs1, Run* mcs2) {

	std::vector<double> dagSteps;
	double dagDown = getFillEnd(mcs1,mcs2) + 10.0; // This will break on non-cleaning runs
	//double acUp = getACUp(mcs1);
	//if (dagDown < acUp) {
	//	dagDown = acUp + 5.0;
		//printf("dagDown = %f, acUp = %f\n",dagDown,acUp);
	//}
	
	//double stepTime = mcs2->getTagBitEvt(1<<0, dagDown+0.25, 1);
	double stepTime = mcs2->getTagBitEvt(1<<0, dagDown-2.0, 1);
	//double tdMove = mcs1->getTagBitEvt((1<<4), stepTime, 1);
	double stop;
	if(stepTime < 0) {
		printf("%f Error! Could not find dagger step. Returning.\n",dagDown);
		dagSteps.clear();
		return dagSteps;
	}

	double lastStep = 0;
	do {
		dagSteps.push_back(stepTime);
		stop = mcs2->getTagBitEvt(1<<0, dagSteps.back() + 0.25, 0);
		stepTime = mcs2->getTagBitEvt(1<<0, stop + 0.25, 1);
		
		/*if (stepTime - dagSteps.back() > 50) { 
			// This is a check for timing 
			// If we have a big dagger step (longer than 50), we want to 
			// double check that there's nothing weird if we project in the 
			// future.
			double stopTmp  = mcs2->getTagBitEvt(1<<0, stepTime + 0.25, 0);
			double startTmp = mcs2->getTagBitEvt(1<<0, stopTmp + 0.25, 1);
			if ((stopTmp == -1.0) || (startTmp == -1.0) || (stopTmp > startTmp)) {
				break;
			}
		}*/
		//if (stepTime - dagSteps.back() > 201) { // Sometimes long holds break things.
		//	break;
		//}
	//------------------------------------------------------------------
	// This while loop tends to mess us up in random things. 
	// For production, want to make sure the last step works.
	// I think the timer was for deep dagger cleaning.
	//------------------------------------------------------------------
	//} while(stepTime != -1.0 && stop != -1.0 && stepTime - dagSteps.back() < 101);
	} while(stepTime != -1.0 && stop != -1.0 && stepTime - dagSteps.back() < 50);
    
    if(stop < 0) {
		printf("%f Error! Could not find dagger step. Returning.\n",dagDown);
        //printf("Error! Could not find dagger step. Returning.\n");
        dagSteps.clear();
		return dagSteps;
    }
	// We've broken out of the loop by getting a "bad" step.
	
	
	/*std::vector<input_t> cts = mcs1->getCounts(
		[](input_t x)->input_t{return x;},
		[dagSteps](input_t x)->bool{return x.realtime > dagSteps.back();}
	);*/
	// Load last timing
	std::vector<input_t> cts = mcs1->getCounts(
		[](input_t x)->input_t{return x;},
		[dagSteps](input_t x)->bool{return x.realtime > dagSteps.back();}
	);
	//dagSteps.push_back(cts.back().realtime);
	// Testing here to get 
	//if(stepTime > 0.0 && stepTime-dagSteps.back() < 50) {
	if(stepTime > 0.0 && (cts.back().realtime - stepTime) < 10) {// && stepTime-dagSteps.back() < 50) {
		dagSteps.push_back(stepTime);
	} else {
		// If there are no steps use last count in MCS1
		//std::vector<input_t> cts = mcs1->getCounts(
		//	[](input_t x)->input_t{return x;},
		//	[dagSteps](input_t x)->bool{return x.realtime > dagSteps.back();}
		//);
		dagSteps.push_back(cts.back().realtime);
	}
	return dagSteps;
}

/* Find times the dagger stops moving  */
std::vector<double> dagStops(Run* mcs1, Run* mcs2, std::vector<double> dips) {

	std::vector<double> dagEnds;
	
	if (dips.size() > 1) { // Need at least 2 movements
		for (size_t ii = 0; ii < dips.size(); ii++) {
			// Ignore the last step since that's just EOR
			double stop = mcs2->getTagBitEvt(1<<0, dips.at(ii) + 0.25, 0);
			dagEnds.push_back(stop);
		}
	} else if (dips.size() == 1) { // assume instantaneous motion for 1-dip
		dagEnds = dips;
	}
	return dagEnds;
		
 }
	

/* Get H-GX hits while the beam (LX current) is on -- Not Online Yet!!!*/
std::vector<double> getLXOn(Run* mcs1, Run* mcs2, EMS* ems0) {
	
	
	std::vector<double>hits;
	// Load all LX counts from LX current
	std::vector<input_e> LXhits = ems0->getCounts(
											[](input_e x)->input_e{return x;},
											[](input_e x)->bool{return (x.dev == 9001) && (x.ch == 1);} 
										);
	// Initialize looping variables
	// lxTime is imported as an int
	double pulseStart = 0.0;
	double lxTime = 0.0;
	input_e event1;
	input_e event2;
	
	double LXval;
	double iniTime = (double)(LXhits[1].time);
	double prevTime;
	double postTime;
	//printf("%d %f\n", LXhits[1].time, (double)LXhits[1].time);
	//printf("%lu\n", LXhits.size());
	// Loop through events from LXhits vector
	//for(auto stepIt = LXhits.begin() + 1; stepIt < LXhits.end(); stepIt++) {
	for (int i = 1; i < LXhits.size()-1; i++) {	
		// Load our two events
		//printf("%d\n",i);
		
		prevTime = (double)(LXhits[i].time)-iniTime;
		postTime = (double)(LXhits[i+1].time)-iniTime;
		LXval = LXhits[i].val;
		//event1 = (double) (*(stepIt - 1).;
		//event2 = *(stepIt);
		//printf("%f,%f,%f\n", prevTime, postTime, LXval);
		// If LX is zero, we don't get events! Arbitrary spacing (10s) necessary here
		// LX value should be ~10. But it should definitely be more than 1.
		if (/*((postTime - prevTime) < 20.0) && */(LXval > 2.0)){
		//if (((event2.time - event1.time) > 10) && (event1.val > 1.0)) {
			
			printf("%f %f\n,", pulseStart, LXval);
			hits.push_back(pulseStart);
			pulseStart = mcs2->getTagBitEvt(1<<1, pulseStart + 1.0, 1);
			if (pulseStart < prevTime) {
				pulseStart = prevTime;
			}
			// Error catch
			if (pulseStart < 0.0) {
				hits.clear();
				break;
			}
		}
	}
	
	// other error catch
	if (hits.front() >= hits.back()) {
		hits.clear();
	}
	
	return hits;//	if pulseSt
		
	//	lxTime = (double)(event2.time);
	//	while(pulseStart < 10.1 && pulseStart >= 0.0) {
	//		hits.push_back(pulseStart);
//			pulseStart = mcs2->getTagBitEvt(1<<1, pulseStart + 2.0, 1);
			//pulseStart = mcs2->getTagBitEvt(1<<1, pulseStart + 2.0, 1);
	//	}
}
	/* ALL H-GX*/
	// Initialize timing with getFillEnd function, and feed
	//double fillEnd = getFillEnd(mcs1,mcs2);
	//PLACEHOLDER
	//double pulseStart = 0.0;
	//double fillEnd = 9999.0;
	
	//std::vector<double> hits;
	
	// loop through all H-GX hits
	//while(pulseStart < fillEnd && pulseStart >= 0.0) {
	//	hits.push_back(pulseStart);
	//	pulseStart = mcs2->getTagBitEvt(1<<1, pulseStart + 2.0, 1);
		//pulseStart = mcs2->getTagBitEvt(1<<1, pulseStart + 2.0, 1);
	//}
//	if(pulseStart < 0.0 || (hits.front() >= hits.back())) {
//		hits.clear();
//	}
	//double fillEnd = beamHits.back() + 3.0 < 300.0 ? beamHits.back() : *(beamHits.end()-2);
		
	//printf("%f,%f\n", hits.front(), hits.back());
//	return hits;
	
	
//}

/* Convert a vector of counts to a "coincidence". 
 * This is the same algorithm as real (telescoping) coincidences, but only triggering
 * on one PMT. NB: this is defined here (rather than in RUN) so that 
 * the coincidence settings can change on a channel to channel basis. */
std::vector<coinc_t> getSelfCoincs(std::vector<input_t> &data, double iniWin, double telWin, int numPH) {
	// Declare initial coincidence settings
	std::vector<coinc_t> coinc;
	coinc_t evt;
	evt.time = 0;
	evt.realtime = 0.0;
	evt.pmt1 = 0;
	evt.pmt2 = 0; // PMT2 will always be zero for a self coinc.
	evt.prompt = 0;
	evt.length = 0;
	if (data.size() < 3) { // Should use .empty() but there's a bug and I don't want to do that...
		return coinc; // need at least 3 total photons to use coincidence algorithm
	}
	
	int peSum = 0;
	int promptSum = 0;
	// looping through the imported data (which should already be in 
	// the wanted range/channel
	for(int ii = 0; ii < data.size()-2; ii++) {
		
		// clear previous buffers 
		peSum = 1; // By definition we start with a photon.
		promptSum = 1;
		input_t prevEvt;
		
		// once we've loaded our dataset, search forwards to find coincidence
		for(int cur = ii+1; cur < data.size()-1; cur++) {
			// check the times of consecutive events.
			// If the two are too far apart, break because it's not a coincidence!
			if(data.at(cur).time - data.at(ii).time > (unsigned long) (iniWin/0.8)) {
				//printf("%lu,%lu,%lu\n",data.at(cur).time, (unsigned long)(iniWin/0.8), data.at(cur).time - data.at(ii).time);
				break;
			} 
			
			peSum +=1;
			promptSum += 1;
			prevEvt = data.at(cur);
			int tailIt; 
			// Now integrate through the tail end
			for (tailIt = cur+1; tailIt < data.size(); tailIt++) {
				// If we're outside of the telescoping window
				if ((data.at(tailIt).time - prevEvt.time) > (unsigned long) (telWin/0.8) 
					&& (data.at(tailIt).time -data.at(ii).time > (unsigned long) (iniWin/0.8))) {
					break;
				}
				peSum +=1;
				if ((data.at(tailIt).realtime - data.at(ii).realtime) < telWin) {
					promptSum += 1;
				}
				prevEvt = data.at(tailIt);
			}
			// Check if we've gone above threshold
			if (promptSum >= numPH) {
			//if (peSum >= numPH) {
				evt.realtime = data.at(ii).realtime;
				evt.time = data.at(ii).time;
				evt.pmt1 = peSum;
				evt.prompt = promptSum;
				evt.length = prevEvt.realtime - data.at(ii).realtime;
				coinc.push_back(evt);
				
				// And move to end of the tail
				ii  = tailIt-1;
				cur = tailIt-1;
			} 
			break;
		}	
	}
	
	return coinc;
}

/* Filter an input_t vector to put in an arbitrary deadtime */
std::vector<input_t> imposeDeadtime(std::vector<input_t> &cts, double deadtime) {
	
	std::vector<input_t>  ctsFilter;
	double prevT = 0.0;
	// cts should already be sorted by realtime.
	for (auto cIt = cts.begin(); cIt < cts.end(); cIt++ ) {
		if (cIt->realtime - prevT < deadtime) {
			continue;
		}
		cIt->deadtime = deadtime;
		ctsFilter.push_back(*cIt);
		prevT = cIt->realtime;
	}	
	return ctsFilter;
	
}


/* Filter a coinc_t vector to remove events with a different profile,
 * I'm calling them "electric noise". These are very (unnaturally) fast. 
 * Filtering should remove:
 *  len <= (nPh-PE)*deadt (which itself is 1.5*dt)*/
std::vector<coinc_t> removeElectricNoise_coinc(std::vector<coinc_t> &coinc,double deadt, int PE) {
	// We care about the _size_ of the filtered vector and the 
	// correction we'll get due to RDE (effective deadtime). 
	// Clone filters
	std::vector<coinc_t> filter;
	if (coinc.size()==0) {return filter;} 
	// Here we're going to
	// Go through each event and see if it fits in our intended window 
	for (auto cIt = coinc.begin(); cIt < coinc.end(); cIt++) {
		if (cIt->length >= (((double)(cIt->pmt1+cIt->pmt2)-PE))*deadt) {
			filter.push_back(*cIt);
		} 
	}
	return filter;
}

std::vector<coinc_t> getElectricNoise_coinc(std::vector<coinc_t> &coinc, double deadt, int PE, double tail) {
	
	std::vector<coinc_t> filter;
	if (coinc.size() == 0) {return filter;}
	for (auto cIt = coinc.begin(); cIt < coinc.end(); cIt++) {
		if ((cIt->length - tail) <= (((double)(cIt->pmt1+cIt->pmt2) - PE))*deadt) {
			filter.push_back(*cIt);
		}
	}
	return filter;
}
	
/* Removing photons that are associated with "electric noise" events */
std::vector<input_t> removeElectricNoise_sing(std::vector<input_t> &cts, std::vector<coinc_t> &coinc, double deadt, int PE) {
	
	std::vector<input_t> filter;
	
	double window = 1000*NANOSECOND; // 1us window around "bad" coincidence events
	
	// Checks that we loaded counts
	if (cts.size() == 0) {return filter;}
	if (coinc.size() == 0) {return cts;} // Don't filter if we loaded counts w/o UCN
	
	// Loop through counts and coincidence events
	double wMin = cts.back().realtime; // Find time windows
	double wMax = cts.back().realtime; // Without any additional info, assume last bin
	size_t cInd = 0;
	// Loop through the coincidences and find the first event outside window
	for (cInd = 0; cInd < coinc.size(); cInd++) {
		// We want to remove "bad" events
		int peSum = coinc.at(cInd).pmt1 + coinc.at(cInd).pmt2;
		if (coinc.at(cInd).length < ((double)(peSum-PE)*deadt)) { // first "bad" event
			wMin = coinc.at(cInd).realtime - window; // Minimum time
			wMax = coinc.at(cInd).realtime + coinc.at(cInd).length + window; // Maximum time
			// An event immediately after might still be "electrical noise".
			if ((cInd < coinc.size() - 2) && (coinc.size() > 1)) { // If we're not on the last coinc...
				size_t start = cInd; // Two indices for scanning
				size_t end   = cInd+1;
				
				while((coinc.at(end).realtime - (coinc.at(start).realtime+coinc.at(start).length)) < window) { // Too close!
					start++;
					end++;
					if (end == coinc.size() - 1) {
						start++;
						break; // If "end" is the last then break and bin out all remaining coincidences.
					}
				}
				wMax = coinc.at(start).realtime + coinc.at(start).length + window; // Maximum time
			}
			break; // Only care about first event
		}
	}

	for (size_t ii = 0; ii < cts.size(); ii++) {
		// First check: is this count "above" the coinc
		input_t temp = cts.at(ii);
		
		// Want to remove single photons inside window
		if (cts.at(ii).realtime < wMin) { // wMin is too low
			filter.push_back(temp);
			continue; // In this case, add a count and go on
		}
		// Else, check if we're in the window
		if (cts.at(ii).realtime < wMax) {
			// Is this the first one in the window?
			if (ii == 0) { // If the first photon is the first in a coincidence
				temp.deadtime = wMax - cts.at(ii).realtime;
				filter.push_back(temp);
				continue;
			}
			if ((cts.at(ii).realtime > wMin) && (cts.at(ii-1).realtime < wMin)) {
				temp.deadtime  = wMax - cts.at(ii).realtime; // put the deadtime as the gap
				filter.push_back(temp); // Put this event in
			}
			continue; // And continue onwards
		}
		// Now we need to find the next index, since we're past wMax.
		// This is the same as above.
		for (size_t cInd2 = cInd+1; cInd2 < coinc.size(); cInd2++) {			
			// We want to remove "bad" events
			int peSum = coinc.at(cInd2).pmt1 + coinc.at(cInd2).pmt2;
			if (coinc.at(cInd2).length < ((double)(peSum-PE)*deadt)) { // first "bad" event
				wMin = coinc.at(cInd2).realtime - window; // Minimum time
				wMax = coinc.at(cInd2).realtime + coinc.at(cInd2).length + window; // Maximum time
				// An event immediately after might still be "electrical noise".
				if ((cInd2 < coinc.size() - 2) && (coinc.size() > 1)) { // If we're not on the last coinc...
					size_t start = cInd2; // Two indices for scanning
					size_t end   = cInd2+1;
					while((coinc.at(end).realtime - (coinc.at(start).realtime+coinc.at(start).length)) < window) { // Too close!
						start++;
						end++;
						if (end == coinc.size()-1) {
							start++;
							break; // If "end" is the last then break and bin out all remaining coincidences.
						}
					}
					wMax = coinc.at(start).realtime + coinc.at(start).length + window; // Maximum time
				}
				cInd = cInd2; // And set to the new coincidence buffer
				break; // Only care about first event
			}
		}
		filter.push_back(temp); // The first photon in a flurry is added
	}
	return filter;
}



/* Filter an input_t vector to put in a deadtime */
std::vector<input_t> removeHighUCN_sing(std::vector<input_t> &cts, std::vector<coinc_t> &coinc, double window, int maxPE) {

	std::vector<input_t> filter;
	
	// Checks that we loaded counts
	if (cts.size() == 0) {return filter;}
	if (coinc.size() == 0) {return cts;} // Don't filter if we loaded counts w/o UCN
	// Loop through counts and coincidence events
	double wMin = cts.back().realtime; // Find time windows
	double wMax = cts.back().realtime; // Without any additional info, assume last bin
	size_t cInd = 0;
	// Loop through the coincidences and find the first event with more than maxPE
	for (cInd = 0; cInd < coinc.size()-1; cInd++) {
		if (coinc.at(cInd).pmt1 + coinc.at(cInd).pmt2 > maxPE) { // First event with too many 
			wMin = coinc.at(cInd).realtime; // Minimum time
			size_t start = cInd; // temp indices
			size_t end   = cInd+1;
			//printf("Start,End = %d,%d,%d\n",start,end,coinc.at(cInd).pmt1+coinc.at(cInd).pmt2);
			while ((coinc.at(end).realtime - (coinc.at(end).realtime+coinc.at(start).length)) < window) { // Too close!
				if (end == coinc.size()-1) {
					break;
				} // Safety first!
				start += 1;
				end   += 1;			
				//printf("Start,End = %d,%d\n",start,end);
			}
			if (start != coinc.size() - 1) {
				wMax = coinc.at(start).realtime+coinc.at(start).length+window;
			}else {
				wMax = cts.back().realtime;
			}
			break;
		}
	}
	//printf("wMin,wMax = %f,%f\n",wMin,wMax);
	int peSum = coinc.at(cInd).pmt1 + coinc.at(cInd).pmt2;
	//printf("%d\n",peSum);
	filter.push_back(cts.at(0)); // Put in the first photon
	for (size_t ii = 1; ii < cts.size(); ii++) {
		// First check: is this count "above" the coinc
		input_t temp = cts.at(ii);
		// Find instantaneous rate:
		int nCts = 0;
		for (size_t jj = ii-1; jj > 0; jj--) {
			//printf("%lu,%lu\n",jj,ii);
			if (cts.at(ii).realtime - cts.at(jj).realtime > RATEWINDOW) {
				break;
			}
			nCts += 1;
		}
		//printf("nCts = %d\n",nCts);
		temp.rate = (double)nCts / RATEWINDOW;
		if (cts.at(ii).realtime < wMin) { // wMin is too low
			filter.push_back(temp);
			continue; // In this case, add a run and go on
		}
		// By definition wMin is higher than us -- are we in the window?
		if (cts.at(ii).realtime < wMax) {
			// Is this the first one in the window?
			if ((cts.at(ii).realtime > wMin) && (cts.at(ii-1).realtime < wMin)) {
				// put the deadtime as the gap we use:
				temp.deadtime  = wMin - wMax;
				filter.push_back(temp);
			}		
			continue; // And continue onwards
		}
		// Now we need to find the next index:
		for (size_t cInd2 = cInd; cInd2 < coinc.size()-1; cInd2++) { // New temp value
			wMin = coinc.at(cInd2).realtime;
			if (coinc.at(cInd2).pmt1 + coinc.at(cInd2).pmt2 > maxPE) { // First event with too many 
				wMin = coinc.at(cInd2).realtime; // Minimum time
				size_t start = cInd2; // temp indices
				size_t end   = cInd2+1;
				while ((coinc.at(end).realtime - (coinc.at(end).realtime+coinc.at(start).length)) < window) { // Too close!
					if (end == coinc.size()-1) {
						break;
					} // Safety first!
					start += 1;
					end   += 1;			
					
				}
				if (start != coinc.size() - 1) {
					wMax = coinc.at(start).realtime+coinc.at(start).length+window;
				} else {
					wMax = cts.back().realtime;
				}
				cInd = cInd2;
				break;
			}
		}
		filter.push_back(temp); // The first photon in a flurry is added
	}
	return filter;
}

/* Filter a coinc_t vector to put in a deadtime */
std::vector<coinc_t> removeHighUCN_coinc(std::vector<coinc_t> &coinc, double window, int maxPE) {
	std::vector<coinc_t> filter;
	if (coinc.size() == 0) {return filter;}
	// Much easier to remove here -- just go through and remove any events
	// with more than 100 PE (should be ~1 per run)
	for (auto cIt = coinc.begin(); cIt < coinc.end(); cIt++) {
		// First check if we're too high
		if ((cIt->pmt1 + cIt->pmt2) > maxPE) {
			// Press forwards
			for (auto dIt = cIt; dIt < coinc.end(); dIt++) {
				// Don't want too many UCN close together
				if ((dIt->realtime) > cIt->realtime+cIt->length+window) {
					cIt = dIt-1;
					break;
				}
			}
			continue;
		}
		filter.push_back(*cIt);
	}
	return filter;
}

/* Filter the output of 2 detectors to remove any cross-talk. 
 * e.g. dagger vs. active cleaner */
std::vector<input_t> antiCoincidence(Run* mcs1, std::vector<input_t> &cts1, std::vector<input_t> &cts2) {
	
	// Declare variables
	std::vector<input_t> aCoincs;
	/*input_t evt;
	evt.time = 0;
	evt.realtime = 0.0;
	evt.ch = 0;
	evt.tag = 0;*/
	// Initialize counts 2 iterator
	//input_t* ctsIt = *cts2.begin();
	int window = mcs1->getPeSumWindow();
	
	
	int ctsInd = 0;
	// Want to "antifilter" cts2 against cts1 	
	for(auto ii = cts1.begin(); ii < cts1.end()-1; ii++) {
		
		// Make sure we're on the right ctsIt
		while ((ii->realtime) > (cts2[ctsInd].realtime - window)) {
			if (ctsInd < cts2.size() - 1) {
				ctsInd++;
			} else {
				break;
			}
		}
		
		if (abs(ii->realtime - cts2[ctsInd].realtime) < window) {
			continue;
		} else {
			aCoincs.push_back(cts2[ctsInd]); 
		}
	}
	return aCoincs;
	
}

/* Find every H-GX pulse that can go into our trap */
std::vector<double> hMinGxHits(Run* mcs1, Run* mcs2, double fillEnd) {
	
	double pulseStart = 0.0;
	std::vector<double> hits;
	
	// loop through all H-GX hits
	while(pulseStart < fillEnd && pulseStart >= 0.0) {
		hits.push_back(pulseStart);
		pulseStart = mcs2->getTagBitEvt(1<<1, pulseStart + 2.0, 1);
	}
	if(pulseStart < 0.0 || (hits.front() >= hits.back())) {
		hits.clear();
	}
	
	return hits;
}

/* Fit the fill of our runs. 
 * The fitFill function loads time constant data from an external source.
 * It generates two output files -- one for MAD and one for Exp. */
void fitFill(Run* mcs1, Run* mcs2, const char* filename, int chan) {
	
	// Load detector channel
	Run* detMCS;
	int detChan;
	if (chan <= 5) {
		detMCS = mcs1;	
		detChan = chan + detMCS->getMCSOff();
	} else if ((5 < chan) && (chan <= 10)) {
		detMCS = mcs2;
		detChan = (chan-5) + detMCS->getMCSOff();
	}
	int runNo = mcs1->getRunNo();
	
	// load time constant data
	FILE* tc_file = fopen(filename, "r");
	char line[72];
	double offT; 
	double kap1;
	double kap2;
	while (fgets(line,72,tc_file)) {
		const char* minRun = strtok(line,",");
		const char* maxRun = strtok(NULL,",");
		if ((std::stoi(minRun) <= runNo) && (std::stoi(maxRun) > runNo)) {
			const char* det = strtok(NULL,",");
			if (std::stoi(det) == chan) {
				offT = std::stod(strtok(NULL,","));
				kap1 = std::stod(strtok(NULL,","));
				kap2 = std::stod(strtok(NULL,","));
				break;
			}
		}
	}
	fclose(tc_file);
		
	// Find the fill we need.
	double fillEnd = getFillEnd(mcs1, mcs2);
	
	// Load H-GX hits
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2,fillEnd);
	if(beamHits.size() == 0 || beamHits.back() < 0) {
		printf("Error! Could not find H-GX pulses out to fillEnd!\n");
		return;
	}
	if(*beamHits.begin() >= fillEnd){ 
		printf("Error! Delayed fill! H-GX pulses start after the end of fill!");
		return;
	}
	
	// Load detector counts
	std::vector<input_t> detCts = detMCS->getCounts(
		[fillEnd](input_t x)->input_t{return x;},
		[detChan, fillEnd](input_t x)->bool{return x.ch == detChan && x.realtime < fillEnd;}
	);
	
	TH1D det("histo", "title", (int)(fillEnd), 0,round(fillEnd));
	std::for_each(detCts.begin(), detCts.end(), [&det](input_t x){det.Fill(x.realtime);});
	
	// Generate a fitting function and fit with ROOT
	ExpFill funcMAD;
	funcMAD.setOffsetT(offT);
	funcMAD.setKappa1(kap1);
	funcMAD.setKappa2(kap2);
		
	for(auto it = beamHits.begin(); it < beamHits.end(); it++) {
		if(*it+offT > fillEnd) {
			beamHits.pop_back();
			continue;
		}
		funcMAD.addOffset(*it);
	}

	
	TF1* fitMAD = new TF1("fit", funcMAD, 0.0, fillEnd, beamHits.size());
	for(int i = 0; i < beamHits.size(); i++) {
		fitMAD->SetParameter(i, 800);
		fitMAD->SetParLimits(i, 0.0, 10000);
	}
	det.Fit(fitMAD, "Q"); // make fitting quiet!
	
	// Save MAD structure
	FILE * outfile;
	outfile = fopen("fillDataFixedTC.csv","a");
	printf("FillData - %d,%d,%d,%f,%f\n", runNo,chan,-1,fitMAD->GetChisquare(), (double)fitMAD->GetNDF());
	fprintf(outfile,"%d,%d,%d,%f,%f\n", runNo,chan,-1,fitMAD->GetChisquare(), (double)fitMAD->GetNDF());
	for(int i = 0; i < beamHits.size(); i++) {
		fprintf(outfile,"%d,%d,%d,%f,%f\n", runNo,chan,i,fitMAD->GetParameter(i), fitMAD->GetParError(i));
	}
	fclose(outfile);	
	
	// Fit summed Fill to get time constants
	FillTime funcExp;
	TF1* fitExp = new TF1("fit", funcExp, 0.0, 300.0, 2);
	// Set Limits for parameters
	fitExp->SetParameter(0, 10000); // Total number of UCN
	fitExp->SetParameter(1, 70.0); // Kappa
	fitExp->SetParLimits(1,0.0, 300.0);
		
	det.Fit(fitExp, "Q"); // make fitting quiet!
	
	// Save MAD structure
	outfile = fopen("fillDataExpFit.csv","a");
	fprintf(outfile,"%d,%d,%f,%f\n", runNo,chan,fitExp->GetParameter(1), fitExp->GetParError(1));
	fclose(outfile);						
	printf("Time Constant for channel %d is %f (+- %f) !\n\n", chan,fitExp->GetParameter(1), fitExp->GetParError(1));
	
	delete fitMAD;
	delete fitExp;
}

/* Fit the fill of our runs. 
 * fitFillFree finds the time constants as well.
 * This year we have H-GX data so we can just find the edges */
void fitFillFree(Run* mcs1, Run* mcs2, int chan) {
	
	// Load H-GX hits and fill start/stop
	double fillStart = 0.0;
	double fillEnd = getFillEnd(mcs1, mcs2);
	printf("fillEnd = %f\n",fillEnd);
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2, fillEnd);
	
	if(beamHits.size() == 0 || beamHits.back() < 0) {
		printf("Error! Could not find H-GX pulses out to fillEnd!\n");
		return;
	}
	if(*beamHits.begin() >= fillEnd){ 
		printf("Error! Delayed fill! H-GX pulses start after the end of fill!");
		return;
	}
	
	std::vector<input_t> detCts;
	
	Run* detMCS;
	int detChan;
	if (chan <= 5) {
		detMCS = mcs1;	
		detChan = chan + detMCS->getMCSOff();
	} else if ((5 < chan) && (chan <= 10)) {
		detMCS = mcs2;
		detChan = (chan-5) + detMCS->getMCSOff();
	}
	
	detCts = detMCS->getCounts(
			[fillEnd](input_t x)->input_t{return x;},
			[detChan, fillEnd](input_t x)->bool{return x.ch == detChan && x.realtime < fillEnd;}
	);
	
	TH1D det("histo", "title", (int)fillEnd, fillStart, fillEnd);
	std::for_each(detCts.begin(), detCts.end(), [&det](input_t x){det.Fill(x.realtime);});
	
	// Define a fitting function (of object ExpFillFree
	ExpFillFree func;		
	for(int i = 0; i < beamHits.size(); i++) {
			
		//make sure we're not picking up a pulse after the GV
		if(beamHits[i] + 1.0 > fillEnd) {
			beamHits.pop_back();
		} else {
			func.addOffset(beamHits[i]);
		}
	}
	
	// Define fitting function parameters.	
	TF1* fit = new TF1("fit", func, fillStart, fillEnd, (int)(beamHits.size())+3);
	
	fit->SetParameter(0, 0.5);
	fit->SetParLimits(0, 0.0, 5.0);
	fit->SetParameter(1, 0.5);
	fit->SetParLimits(1, 0.0, 10.0);
	fit->SetParameter(2, 20.0);
	fit->SetParLimits(2, 0.0, 50.0);
	for(int i = 0; i < (int)(beamHits.size()); i++) {
		fit->SetParameter(i+3, 250);
		fit->SetParLimits(i+3, 0.0, 1000.0);
	}
		
	det.Fit(fit,"Q"); // make fitting quiet!
	
	// Output to terminal and file
		
	int runNo = mcs1->getRunNo();
	FILE * outfile;
	outfile = fopen("fillDataFree.csv","a");
		
	printf("FillData - %d,%d,%d,%f,%f\n", runNo, chan,-1, fit->GetChisquare(), (double)fit->GetNDF());
	fprintf(outfile,"%d,%d,%d,%f,%f\n", runNo, chan,-1, fit->GetChisquare(), (double)fit->GetNDF());
	for(int i = 0; i < (int)(beamHits.size()) + 3; i++) {
		printf("FillData - %d,%d,%d,%f,%f\n", runNo,chan, i, fit->GetParameter(i), fit->GetParError(i));
		fprintf(outfile,"%d,%d,%d,%f,%f\n", runNo,chan, i, fit->GetParameter(i), fit->GetParError(i));
	}
	fclose(outfile);	
	
	delete fit;
}

/* Fit the fill of our runs, in multiple detectors.
   This year we have H-GX data so we can just find the edges */
void fitFillMulti(Run* mcs1, Run* mcs2) {
		
	//double fillEnd = 300.0; // fill with Roundhouse
	double fillEnd = getFillEnd(mcs1,mcs2);
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2, fillEnd);
	if(beamHits.size() == 0 || beamHits.back() < 0) {
		printf("Error! Could not find H-GX pulses out to fillEnd!\n");
		return;
	}
	if(*beamHits.begin() >= fillEnd){ 
		printf("Error! Delayed fill! H-GX pulses start after the end of fill!");
		return;
	}
	
	// Load all possible filling parameters
	std::vector<input_t> d1Cts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{return x;},
		[fillEnd](input_t x)->bool{return x.ch == 1 && x.realtime < fillEnd;}
	);
	std::vector<input_t> d2Cts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{return x;},
		[fillEnd](input_t x)->bool{return x.ch == 2 && x.realtime < fillEnd;}
	);
	std::vector<input_t> olCts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{return x;},
		[fillEnd](input_t x)->bool{return x.ch == 3 && x.realtime < fillEnd;}
	);
	std::vector<input_t> baCts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{return x;},
		[fillEnd](input_t x)->bool{return x.ch == 4 && x.realtime < fillEnd;}
	);
	std::vector<input_t> spCts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{return x;},
		[fillEnd](input_t x)->bool{return x.ch == 5 && x.realtime < fillEnd;}
	);
	
	// initialize ROOT histograms and count them
	TH1D d1("histo", "title", (int)(fillEnd), 0, fillEnd);
	TH1D d2("histo", "title", (int)(fillEnd), 0, fillEnd);
	TH1D ol("histo", "title", (int)(fillEnd), 0, fillEnd);
	TH1D ba("histo", "title", (int)(fillEnd), 0, fillEnd);
	TH1D sp("histo", "title", (int)(fillEnd), 0, fillEnd);
	std::for_each(d1Cts.begin(), d1Cts.end(), [&d1](input_t x){d1.Fill(x.realtime);});
	std::for_each(d2Cts.begin(), d2Cts.end(), [&d2](input_t x){d2.Fill(x.realtime);});
	std::for_each(olCts.begin(), olCts.end(), [&ol](input_t x){ol.Fill(x.realtime);});
	std::for_each(baCts.begin(), baCts.end(), [&ba](input_t x){ba.Fill(x.realtime);});
	std::for_each(spCts.begin(), spCts.end(), [&sp](input_t x){sp.Fill(x.realtime);});
	
	// Assuming idealized fills. 
	ExpFillFree func;
	
	// Figure out how many counts we have.
	// If the last H-GX hit is after the end of fill, remove from fill fitting.
	for(auto it = beamHits.begin(); it < beamHits.end(); it++) {
		if(*it+3.0 > fillEnd) {
			beamHits.pop_back();
			continue;
		}
		func.addOffset(*it);
	}
	int nHits = beamHits.size();

	int i;
	int runNo = mcs1->getRunNo();
	
	TF1* fit = new TF1("fit", func, 0.0, fillEnd, beamHits.size());
	for(i = 0; i < beamHits.size(); i++) {
		fit->SetParameter(i, 800);
		fit->SetParLimits(i, 0.0, 10000);
	}
	sp.Fit(fit, "Q"); // make fitting quiet!
	FILE * outfile;
	outfile = fopen("fillDataMulti.csv","a");
	printf("FillData - %d,%d,%f,%f\n", runNo, -1, fit->GetChisquare(), (double)fit->GetNDF());
	fprintf(outfile,"%d,%d,%f,%f\n", runNo, -1, fit->GetChisquare(), (double)fit->GetNDF());
	for(i = 0; i < beamHits.size(); i++) {
	//	printf("FillData - %d,%d,%f,%f\n", runNo, i, fit->GetParameter(i), fit->GetParError(i));
		fprintf(outfile,"%d,%d,%f,%f\n", runNo, i, fit->GetParameter(i), fit->GetParError(i));
	}
	fclose(outfile);	
	char fName[256];
	sprintf(fName, "summaryPlots/FitFill%05d.root", runNo);
	sp.SaveAs(fName);
	delete fit;
}

/* Fit Background to Multiexponential */
void fitBkgUnload(Run* mcs1, Run* mcs2, EMS* ems) {
	
	// Load tagbit data (dagger movement)
	int runNo = mcs1->getRunNo();
	int cMode = mcs1->getCoincMode();
	
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
	//std::vector<double> dagStop  = dagStops(mcs1,mcs2,dagSteps);
	
	// Make sure the run actually loaded...
	if(dagSteps.empty() || (dagSteps.end() <= dagSteps.begin())) { 
		fprintf(stderr, "Skipping run %d for no Dagger Movement\n", runNo);
		return;
	}

	// Find background times
	// During production running, the last 50s should have TD open, dagger down bkgs
	double countStart;
	double countEnd;
	double bkgStart;
	double bkgEnd = *(dagSteps.end()-1);
	if (dagSteps.size() > 2) { 
		countStart = mcs2->getTagBitEvt(1<<0,*(dagSteps.end()-2)+0.25, 0) + 1.0; // Dagger stops moving, add an extra second
		countEnd = mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);
		bkgStart = mcs1->getTagBitEvt(1<<4, *(dagSteps.end()-2), 1) + 5.0;// 1<<4 is Cat Door, add extra 5 seconds to be safe
	} else {
		countStart = *(dagSteps.begin()-1) - 190; // 150s unload on last dip, assume 10s for movement.
		countEnd = *(dagSteps.end()-1) - 50; // 50
		bkgStart = *(dagSteps.end()-1) - 40;
		fprintf(stderr, "Using FIXED TIMINGS on run %d for background!!!!!!!\n", runNo);
	}
	if (countStart < 0. || countStart >= countEnd) { 
		fprintf(stderr,"Skipping run %d for counting timing glitch\n",runNo);
		return;
	}
	if (bkgStart < 5.0 || bkgStart >= bkgEnd) {
		fprintf(stderr, "Skipping run %d for background timing glitch\n", runNo);
		return;
	}
	//printf("%f,%f,%f,%f\n",countStart,countEnd,bkgStart,bkgEnd);
	// Amount of time to exclude due to trap door movement
	double tdSlop = bkgStart - countEnd;
	
	//------------------------------------------------------------------
	// Could potentially optimize by loading a TC file here. I'm lazy.
	//------------------------------------------------------------------
	
	// Determine which dagger pair we want to use
	Run* dagMCS;
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	// Load background (End of Run) counts 
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
	// Now let's get background rate guesses:
	double bkg1 = ((double) bkgDag1.size()) / (bkgEnd-bkgStart);
	double bkg2 = ((double) bkgDag2.size()) / (bkgEnd-bkgStart);
	double bkgC = ((double) bkgDagC.size()) / (bkgEnd-bkgStart);
		
	// Load other counts vectors.
	// Note that there's a little trickyness going on to exclude TD movement
	std::vector<input_t> ctsDag1 = dagMCS->getCounts(
		[countStart,bkgStart,tdSlop](input_t x)->input_t{
										input_t y = x;
										y.realtime -= countStart;
										//y.realtime < bkgStart ? y.realtime -= countStart : y.realtime -= (countStart+tdSlop);
										//if (x.realtime < bkgStart) {y.realtime -= countStart;}
										//else {y.realtime -= (countStart+tdSlop);}
										return y;
									},
		[c1,countStart,countEnd,bkgStart,bkgEnd](input_t x)->bool{
													return (x.ch == c1 && 
															((x.realtime > countStart && x.realtime < countEnd) || 
															(x.realtime > bkgStart && x.realtime < bkgEnd)));
												}
	);
	std::vector<input_t> ctsDag2 = dagMCS->getCounts(
		[countStart,bkgStart,tdSlop](input_t x)->input_t{
										input_t y = x; 
										y.realtime -= countStart;
										//y.realtime < bkgStart ? y.realtime -= countStart : y.realtime -= (countStart+tdSlop);
										//if (x.realtime < bkgStart) {y.realtime -= countStart;}
										//else {y.realtime -= (countStart+tdSlop);}
										return y;
									},
		[c2,countStart,countEnd,bkgStart,bkgEnd](input_t x)->bool{
													return (x.ch == c2 && 
															((x.realtime > countStart && x.realtime < countEnd) || 
															(x.realtime > bkgStart && x.realtime < bkgEnd)));
												}
	);
	std::vector<coinc_t> ctsDagC = dagMCS->getCoincCounts(
		[countStart,bkgStart,tdSlop](coinc_t x)->coinc_t{
										coinc_t y = x; 
										y.realtime -= countStart;
										//y.realtime < bkgStart ? y.realtime -= countStart : y.realtime -= (countStart+tdSlop);
										//if (x.realtime < bkgStart) {y.realtime -= countStart;}
										//else {y.realtime -= (countStart+tdSlop);}
										return y;
									},
		[c1,countStart,countEnd,bkgStart,bkgEnd](coinc_t x)->bool{
													//return ((x.realtime > countStart && x.realtime < countEnd) || 
													//		(x.realtime > bkgStart && x.realtime < bkgEnd));
													return (x.realtime > countStart && x.realtime < bkgEnd);
												}
	);
	
	// We're doing root fitting. 
	// Let's start by filling histograms
	// Auto-set the number of bins, based on the amount of time in the unload (excluding TD movement)
	//double bins = bkgEnd-countStart-tdSlop;
	double bins = bkgEnd-countStart;
	// Need to rescale singles to treat as "Poisson"
	//double scale1 = sqrt(((double) ctsDagC.size()) /((double) ctsDag1.size()));
	//double scale2 = sqrt(((double) ctsDagC.size()) /((double) ctsDag2.size()));
	double scale1 = 1.;
	double scale2 = 1.;
	TH1D histC("histoC", "title", (int)(bins),0,(int)(bins));
	std::for_each(ctsDagC.begin(),ctsDagC.end(), [&histC](coinc_t x){histC.Fill(x.realtime);});
	TH1D hist1("histo1", "title", (int)(bins),0,(int)(bins));
	std::for_each(ctsDag1.begin(),ctsDag1.end(), [&hist1,scale1](input_t x){hist1.Fill(x.realtime,scale1);});
	TH1D hist2("histo2", "title", (int)(bins),0,(int)(bins));
	std::for_each(ctsDag2.begin(),ctsDag2.end(), [&hist2,scale2](input_t x){hist2.Fill(x.realtime,scale2);});
	
	
	
	// Now we load our functional object:
	size_t nTC = 3; // For now, 3-exponential
	std::vector<int> enTC;
	for (int ii = 0; ii < nTC; ii++) { enTC.push_back(pow(10,ii));}  
	
	hist2.SaveAs("test.root");
	// For Coincidences first...
	MultiExpBkg funcC(nTC); 
	funcC.setBackground(bkgC);
	double peakC = histC.GetMaximum();
	// A little bit of hardcoding here, unfortunately.
	TF1* fitC = new TF1("fitC",funcC,0,bins,2*nTC+1);
	fitC->SetParameter(0,bkgC);
	fitC->SetParLimits(0,0.,100.);
	
	//fitC->SetParameter(1,peakC);
	//fitC->SetParameter(3,peakC/10);
	//fitC->SetParameter(5,peakC/100);
	// If we had a loaded file, there would be a good way to avoid this.
	for (int ii = 0; ii < nTC; ii++) {
		fitC->SetParameter(2*ii+1,peakC/enTC[ii]);
		fitC->SetParLimits(2*ii+1,0.,5*peakC/enTC[ii]);
		fitC->SetParLimits(2*(ii+1),0.,40.);
	}
	// Hard code time constants (guessing these here.)
	fitC->SetParameter(2,8.);
	fitC->SetParameter(4,3.); 
	fitC->SetParameter(6,11.);
	
	
	histC.Fit(fitC,"QL"); // Actual fitting here!
	//printf("%f\n",peak);
	
	// Now do the same thing for Singles. 
	// BUT now we fix to the coincidence time constants
	MultiExpBkg func1(nTC);
	func1.setBackground(bkg1);
	MultiExpBkg func2(nTC);
	func2.setBackground(bkg2);
	double peak1 = hist1.GetMaximum();
	double peak2 = hist2.GetMaximum();
	
	TF1* fit1 = new TF1("fit1",func1,0,bins,2*nTC+1);
	fit1->SetParameter(0,bkg1);
	fit1->SetParLimits(0,0,500);
	TF1* fit2 = new TF1("fit1",func2,0,bins,2*nTC+1);
	fit2->SetParameter(0,bkg2);
	fit2->SetParLimits(0,0,500);
	// If we had a loaded file, there would be a good way to avoid this.
	for (size_t ii = 0; ii < nTC; ii++) {
		fit1->SetParameter(2*ii+1,peak1/enTC[ii]); // Assume 10x coincidence
		fit1->SetParLimits(2*ii+1,0.,5*peak1/enTC[ii]);
		fit1->FixParameter(2*(ii+1),fitC->GetParameter(2*(ii+1))); // Fixing time constants to Coinc.
		fit2->SetParameter(2*ii+1,peak2/enTC[ii]);
		fit2->SetParLimits(2*ii+1,0.,5*peak2/enTC[ii]);
		fit2->FixParameter(2*(ii+1),fitC->GetParameter(2*(ii+1))); // Fixing time constants to Coinc.
	}
	
	hist1.Fit(fit1, "QW"); // Actual fitting here!
	hist2.Fit(fit2, "QW"); 
	
	// Now we save backgrounds (And I guess Time Constants too)!
	FILE* outfile;
	outfile = fopen("multiExpBkg.csv","a");
	printf("MultiExpData - Run %d\n",runNo);
	printf("   Time Constants: \n");
	printf("      1: %f +/- %f\n",fitC->GetParameter(2),fitC->GetParError(2));
	printf("      2: %f +/- %f\n",fitC->GetParameter(4),fitC->GetParError(4));
	printf("      3: %f +/- %f\n",fitC->GetParameter(6),fitC->GetParError(6));
	printf("   Coincidence: Chi2(/NDF) %f (%f)\n",fitC->GetChisquare(),fitC->GetChisquare()/(double)fitC->GetNDF());
	printf("      Background: %f +/- %f\n",fitC->GetParameter(0),fitC->GetParError(0));
	printf("      Compare: %f\n",bkgC);
	printf("      1: %f +/- %f\n",fitC->GetParameter(1),fitC->GetParError(1));
	printf("      2: %f +/- %f\n",fitC->GetParameter(3),fitC->GetParError(3));
	printf("      3: %f +/- %f\n",fitC->GetParameter(5),fitC->GetParError(5));
	printf("   PMT 1: Chi2(/NDF) %f (%f)\n",fit1->GetChisquare(),fit1->GetChisquare()/(double)fit1->GetNDF());
	printf("      Background: %f +/- %f\n",fit1->GetParameter(0),fit1->GetParError(0));
	printf("      Compare: %f\n",bkg1);
	printf("      1: %f +/- %f\n",fit1->GetParameter(1),fit1->GetParError(1));
	printf("      2: %f +/- %f\n",fit1->GetParameter(3),fit1->GetParError(3));
	printf("      3: %f +/- %f\n",fit1->GetParameter(5),fit1->GetParError(5));
	printf("   PMT 2: Chi2(/NDF) %f (%f)\n",fit2->GetChisquare(),fit2->GetChisquare()/(double)fit2->GetNDF());
	printf("      Background: %f +/- %f\n",fit2->GetParameter(0),fit2->GetParError(0));
	printf("      Compare: %f\n",bkg2);
	printf("      1: %f +/- %f\n",fit2->GetParameter(1),fit2->GetParError(1));
	printf("      2: %f +/- %f\n",fit2->GetParameter(3),fit2->GetParError(3));
	printf("      3: %f +/- %f\n",fit2->GetParameter(5),fit2->GetParError(5));
	
	fprintf(outfile,"%d,",runNo);
	fprintf(outfile,"%f,%f,%f,%f,", fitC->GetChisquare(),(double)fitC->GetNDF(),fitC->GetParameter(0),fitC->GetParError(0));
	fprintf(outfile,"%f,%f,%f,%f,", fit1->GetChisquare(),(double)fit1->GetNDF(),fit1->GetParameter(0),fit1->GetParError(0));
	fprintf(outfile,"%f,%f,%f,%f\n",fit2->GetChisquare(),(double)fit2->GetNDF(),fit2->GetParameter(0),fit2->GetParError(0));
	fclose(outfile);
	
	// Since I'm fitting stuff, might as well save time constants for later.
	FILE* tcOutfile;
	tcOutfile = fopen("tcOutputs.csv","a");
	fprintf(tcOutfile,"%d,",runNo);
	for (size_t ii = 0; ii < nTC; ii++) {
		fprintf(tcOutfile,"%f,%f,%f,%f,", fitC->GetParameter(2*ii+1),fitC->GetParError(2*ii+1),fitC->GetParameter(2*(ii+1)),fitC->GetParError(2*(ii+1)));
		fprintf(tcOutfile,"%f,%f,%f,%f,", fit1->GetParameter(2*ii+1),fit1->GetParError(2*ii+1),fit1->GetParameter(2*(ii+1)),fit1->GetParError(2*(ii+1)));
		fprintf(tcOutfile,"%f,%f,%f,%f\n", fit2->GetParameter(2*ii+1),fit2->GetParError(2*ii+1),fit2->GetParameter(2*(ii+1)),fit2->GetParError(2*(ii+1)));
	}
	fclose(tcOutfile);
	// Clear Memory
	delete fitC;
	delete fit1;
	delete fit2;
	
}


/* Fit the last unload and generate an output for ZnS simulations */
void fitLastDip(Run* mcs1, Run* mcs2, TH1D* dipFit) {
	
	// Load coincidence mode for high vs low threshold
	int cMode = mcs1->getCoincMode();
	int runNo = mcs1->getRunNo();
	
	// Determine which dagger run we want to use
	Run* dagMCS;
	if ((cMode % 2) == 1) { 
		dagMCS = mcs1; // Low Thresh.
	} else {
		dagMCS = mcs2; // High Thresh.
	}

	// And get channels
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	// Load dagSteps
	std::vector<double> dagSteps = dagDips(mcs1, mcs2);
	// make sure we're not just doing a 1 step (or aborted run)
	if(dagSteps.size() < 2) {
		return;
	}
	
	// Find timing numbers for last dip only
	double startTime  = *(dagSteps.end() - 2);
	double bkgStart= mcs1->getTagBitEvt(1<<3, *(dagSteps.end()-2), 1);;
	double bkgEnd = *(dagSteps.end()-1);
	
	// For now separate short run.
	if (startTime > 650 || bkgStart > 850) {
		printf("Skipping Long Run!\n");
		return;
	}
	
	// Load coincidence timing counter
	std::vector<coinc_t> coinc;
	
	coinc = dagMCS->getCoincCounts(
		[startTime, bkgStart](coinc_t x)->coinc_t{x.realtime -= startTime; return x;},
		[startTime, bkgStart](coinc_t x)->bool{return (x.realtime > startTime && x.realtime < bkgStart);}
	);
	
	// Make sure we get actual counts
	if(coinc.size() < 20) { 
		return;
	}

	// Initialize, fill, and normalize root histogram
	TH1D hist("histo", "title", 300, 0, 150);
	std::for_each(coinc.begin(), coinc.end(), [&hist](coinc_t x){hist.Fill(x.realtime);});	
	hist.Scale(1.0/coinc.size());
	
	// Load our draining time function
	DrainTime func;
	TF1* fit = new TF1("fit", func, 0.0, 50.0, 3);
	// Set Limits for parameters
	fit->SetParameter(0, 10); // T0
	fit->SetParLimits(0, 0.0, 2000);
	fit->SetParameter(1, 5); // Kappa
	fit->SetParLimits(1,0.0, 10);
		
	hist.Fit(fit, "Q"); // make fitting quiet!
	
	// Output file
	FILE * outfile;
	outfile = fopen("lastDipFit.csv","a");
	printf("FillData - %d,%d,%f,%f\n", runNo, 0, fit->GetParameter(0), fit->GetParError(0));
	fprintf(outfile,"%d,%d,%f,%f\n", runNo, 0, fit->GetParameter(0), fit->GetParError(0));
	printf("FillData - %d,%d,%f,%f\n", runNo, 1, fit->GetParameter(1), fit->GetParError(1));
	fprintf(outfile,"%d,%d,%f,%f\n", runNo, 1, fit->GetParameter(1), fit->GetParError(1));
	fclose(outfile);
	
	// rescale hist and add to overall sum
	hist.Scale(coinc.size());
	dipFit->Add(&hist);

	// save memory by deleting root stuff
	delete fit;
}

//----------------------------------------------------------------------
// These functions compile, but have been absorbed into other things.
// TODO: Update these
//----------------------------------------------------------------------
/* Measure efficiency of the dagger in peak 1*/
void effMeas(Run* mcs1, Run* mcs2) {
	
	struct coinc {
		double t;
		int ch;
	};
	// Find when the first dip happens in mcs1
	double firstDip = mcs1->getTagBitEvt(1<<9, 175, 0);
	firstDip += 40.0;
	std::vector<coinc_t> dCts = mcs1->getCoincCounts( // Should do something about high threshold.
		[firstDip](coinc_t x)->coinc_t{x.realtime -= firstDip; return x;},
		[firstDip](coinc_t x)->bool{return x.realtime > firstDip && x.realtime < firstDip + 20;}
	);
	std::vector<coinc_t> acCts = mcs2->getCoincCounts(
		[firstDip](coinc_t x)->coinc_t{x.realtime -= firstDip; return x;},
		[firstDip](coinc_t x)->bool{return x.realtime > firstDip && x.realtime < firstDip + 20;}
	);
	// Make sure there's no empty coincidences
	std::vector<coinc> allCts;
	if(dCts.empty() || acCts.empty()) {
		printf("Empty! dCts: %ld aCts: %ld\n", dCts.size(), acCts.size());
		return;
	}
	// load and sort our events
	std::for_each(dCts.begin(), dCts.end(), [&allCts](coinc_t x){coinc y = {x.realtime, 1}; allCts.push_back(y);});
	std::for_each(acCts.begin(), acCts.end(), [&allCts](coinc_t x){coinc y = {x.realtime, 2}; allCts.push_back(y);});
	std::sort(allCts.begin(), allCts.end(), [](coinc x, coinc y)->bool{return (x.t < y.t);});
	
	// scan across all events to find the number of coincidences
	int numCoinc = 0;
	for(auto ii = allCts.begin(); ii < allCts.end(); ii++) {
		for(auto jj = ii+1; jj < allCts.end(); jj++) {
			if(jj->t - ii->t > 60000 * NANOSECOND) {
				break;
			}
			if(jj->ch != ii->ch) {
				numCoinc += 1;
				ii = jj;
				break;
			}
		}
	}
	
	printf("Dagger Efficiency: %f - coinc/ac: %d / %lu dag: %lu; predicted %f\n", ((double)numCoinc)/((double)acCts.size()), numCoinc, acCts.size(), dCts.size(), dCts.size() * (double)acCts.size() / (double)numCoinc);
}

/* Estimate the number of UCN in the trap -- presently obsolete becasue
 * I do this in Python. */
void estimateN0(Run* mcs1, Run* mcs2) {
	double fillEnd = getFillEnd(mcs1, mcs2);
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2, fillEnd);
	if(beamHits.empty()) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", mcs1->getRunNo());
		return;
	}
	//double fillEnd = getFillEnd(mcs1,mcs2);
	std::vector<input_t> spCts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 5 && x.realtime < fillEnd;}
	);
	std::vector<input_t> bareCts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 4 && x.realtime < fillEnd;}
	);
	measurement weightSP = expWeightMonVect(spCts, 70.0);
	measurement weightBare = expWeightMonVect(bareCts,70.0);
	
	double start = mcs1->getTagBitEvt(1<<9, 765, 0);
	if(fabs(start - 770) > 2.0) {
		printf("No step at 770!!!\n");
		return;
	}
	double stop = start + 100;
	
	std::vector<coinc_t> dCts = mcs1->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[start, stop](coinc_t x)->bool{return x.realtime > start && x.realtime < stop;}
	);

	printf("Data - %d,%f,%lu,%f,%f,%f,%f\n", mcs1->getRunNo(), start, dCts.size(), weightSP.val, weightSP.err, weightBare.val, weightBare.err);	
	//fitFill(mcs1, mcs2);
}

/* Scan the number of counts in the beginning vs end of run
 * 
 * I think this was for a specific systematic run that we haven't done
 * in 2017 or 2018. */
void fillingScanMeas(Run* mcs1) {
	// Get coincidence coints from Run object
	std::vector<coinc_t> dCts = mcs1->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[](coinc_t x)->bool{return true;}
	);
	// Create a vector of just coincidence hits
	std::vector<unsigned int> cHist;
	for(auto ii = dCts.begin(); ii < dCts.end(); ii++) {
		while(cHist.size() <= (size_t)floor(ii->realtime)) {
			cHist.push_back(0);
		}
		cHist[(size_t)floor(ii->realtime)] += 1;
	}
	// Break apart number of coincidences
	if(cHist.size() > 400) { // First 400 UCN hits
		auto maxBin = std::max_element(cHist.begin(), cHist.begin()+400);
		if(maxBin-10 >= cHist.begin() && maxBin+100 < cHist.end()) {
			unsigned long integral = std::accumulate(maxBin-10, maxBin+100, 0, [](unsigned long sum, unsigned int x){return sum+x;});
			printf("1st: %lu\n", integral);
		}
	}
	if(cHist.size() > 800) { // Remaining UCN hits
		auto maxBin = std::max_element(cHist.begin()+400, cHist.end());
		if(maxBin-10 >= cHist.begin() && maxBin+100 < cHist.end()) {
			unsigned long integral = std::accumulate(maxBin-10, maxBin+100, 0, [](unsigned long sum, unsigned int x){return sum+x;});
			printf("2nd: %lu\n", integral);
		}
	}
}



/* Normalized neutrons for Singles analysis
void extractObservables(Run* mcs1, Run* mcs2, const char* filename) {
	
	// Coinc By Dip+bkg (with error = sqrt coinc)
	// Sing By Dip+bkg (with error = sqrt sing/8)
	
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
	
	// Determine which dagger pair we want to use
	if ((mcs1->getCoincMode() == 1) || (mcs1->getCoincMode() == 2)) {
		dagMCS = mcs1;
	} else if ((mcs1->getCoincMode() == 3) || (mcs1->getCoincMode() == 4)) {
		dagMCS = mcs2;
	}
	int c1 = dagMCS->getCoinC1() + dagMCS->getMCSOff();
	int c2 = dagMCS->getCoinC2() + dagMCS->getMCSOff();
	
	// Load background counts
	std::vector<input_t> bkgDag1;
	std::vector<input_t> bkgDag2;
	std::vector<input_t> bkgDagC;
	bkgDag1 = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c1,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c1 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	bkgDag2 = dagMCS->getCounts(
		[](input_t x)->input_t{return x;},
		[c2,bkgStart,bkgEnd](input_t x)->bool{return (x.ch == c2 && x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	bkgDagC = dagMCS->getCoincCounts(
		[](input_t x)->input_t{return x;},
		[bkgStart, bkgEnd](input_t x)->bool{return (x.realtime > bkgStart && x.realtime < bkgEnd);}
	);			
			
	// Output 2 .csv files that we will use in NormAnalyzer-Coinc.py
	// Start by outputting monitor and backgrounds
	FILE* outfile;
	outfile = fopen("ExtractedObservables.csv", "a");
	printf("Outputting monitor/background data to Ecsv!\n");
	printf("\nMonitor Data -----Run: %d\nfillEnd: %f\nholdEnd: %f\bkgDag1: %lu\nbkgDag2: %lu\nbkgDagC: %lu\nbkgStart: %f\nbkgEnd: %f\n",
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
	FILE* outfileDag;
	outfileDag = fopen("normNByDipSing-dag.csv","a");
	
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
		
		if (startTime >= endTime) {
			printf("Stopping due to timing bug!\n");
			return;
		}
		
		// Divide each step into 10 second slices
		double tempSta = startTime;
		double tempEnd = startTime + sliceTime;
		for (int jj = 0; jj < (int)((endTime - startTime) / sliceTime) + 1; jj++) {
			
			// Make sure we stop at the right step
			if (endTime > tempEnd) {
				tempEnd = endTime;
			}
			if ((tempSta + 0.2) > tempEnd) { // 0.2s timing slop
				break;
			}
			
			// Load counts (singles)
			std::vector<input_t> cts1 = dagMCS->getCounts(
				[](input_t x)->input_t{return x;},
				[c1,tempSta,tempEnd](input_t x)->bool{return (x.ch == c1 && x.realtime > tempSta && x.realtime < tempEnd);}
			);
			std::vector<input_t> cts2 = dagMCS->getCounts(
				[](input_t x)->input_t{return x;},
				[c2,tempSta,tempEnd](input_t x)->bool{return (x.ch == c2 && x.realtime > tempSta && x.realtime < tempEnd);}
			);
			
			// Deadtime correction for singles
			// Step Mean calculation
			double stepMean1 = cts1.size() > 0 ? std::accumulate(cts1.begin(), cts1.end(), 0.0, [](double m, input_t x)->double{return m + x.realtime;}) / (double)cts1.size() : 0.0;
			double stepMean2 = cts2.size() > 0 ? std::accumulate(cts2.begin(), cts2.end(), 0.0, [](double m, input_t x)->double{return m + x.realtime;}) / (double)cts2.size() : 0.0;
						
			// Deadtime Correction
			double time;
			if (cts1.size() != 0) {
				time = floor(cts1.front().realtime);
			} else if (cts2.size() != 0) {
				time = floor(cts2.front().realtime);
			} else {
				printf("No counts loaded! Continuing to next step!\n");
				continue;
			}
			
			// Hardcoded in 20ns with instantaneous rates every s.
			double counts = 0.0;
			double corr1 = 0.0;
			for(auto cIt = cts1.begin(); cIt < cts1.end(); cIt++) {
				if(floor(cIt->realtime) != time) {
					corr1 += (counts/(1.0-counts*20*NANOSECOND) - counts); 
					counts = 0.0;
					time = floor(cIt->realtime);
				}
				counts += 1.0;
			}
			
			counts = 0.0;
			double corr2 = 0.0;
			for(auto cIt = cts2.begin(); cIt < cts2.end(); cIt++) {
				if(floor(cIt->realtime) != time) {
					corr2 += (counts/(1.0-counts*20*NANOSECOND) - counts); // Deadtime correction
					counts = 0.0;
					time = floor(cIt->realtime);
				}
				counts += 1.0;
			}
			
			// Write out dagger data
			printf("\nData -----  Run: %d\numSteps: %d\ntempSta: %f\ntempEnd: %f\nstepMean1: %f\ncts1: %lu\ncorr1: %f\nstepMean2: %f\ncts2: %lu\ncorr2:%f\n",
				   runNo, numSteps, tempSta, tempEnd, stepMean1, cts1.size(), corr1, stepMean2, cts2.size(), corr2);
			fprintf(outfileDag,"%d,%d,%f,%f,%f,%lu,%f,%f,%lu,%f\n",
					runNo, numSteps, tempSta, tempEnd, stepMean1, cts1.size(), corr1, stepMean2, cts2.size(), corr2);
			
			// Increment start and end of step
			tempSta += sliceTime;
			tempEnd += sliceTime;
		}
		numSteps += 1;		
	}
	fclose(outfileDag);	
}*/


/* Output our coincidence root files (in the summaryPlots folder) */
/*void writeCoincHist(Run* run) {
	int runNo = run->getRunNo();
	char fName[256];
	sprintf(fName, "summaryPlots/coincs%05d.root", runNo);
	TH1D dagHist("coinc", "title", 1750, 0, 1750);
	std::vector<input_t> cts = run->getCoincCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return true;});
	for(auto it = cts.begin(); it < cts.end(); it++) {
		dagHist.Fill(it->realtime);
	}
	dagHist.SaveAs(fName);

	sprintf(fName, "summaryPlots/mon%05d.root", runNo);
	TH1D monHist("coinc", "title", 1750, 0, 1750);
	std::vector<input_t> monCts = run->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch==3;});
	for(auto it = monCts.begin(); it < monCts.end(); it++) {
		monHist.Fill(it->realtime);
	}
	monHist.SaveAs(fName);
}*/

/* Normalized neutrons for Singles analysis */
/*void normNByDipSing(Run* mcs1, Run* mcs2) {
	 THIS IS THE OLD VERSION -- USE THIS BECAUSE I SUCK AT VERSION CONTROL
	// Load tagbits: TD and H-GX
	double tdTime = getFillEnd(mcs1, mcs2);;	
	double fillEnd = tdTime;
	std::vector<double> beamHits = hMinGxHits(mcs1, mcs2, fillEnd);
	
	// Make sure the run actually loaded.
	if(beamHits.empty() || (beamHits.end() <= beamHits.begin())) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", mcs1->getRunNo());
		return;
	}
	
	// Figure out when the beam shuts off. 
	// Create vectors of hits in the SP, OLD, and BARE monitors
	//double fillEnd = beamHits.back() + 3.0 < 300.0 ? beamHits.back() : *(beamHits.end()-2);
	std::vector<input_t> oldCts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 3 && x.realtime < fillEnd;}
	);
	std::vector<input_t> bareCts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 4 && x.realtime < fillEnd;}
	);
	std::vector<input_t> spCts = mcs1->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 5 && x.realtime < fillEnd;}
	);
		
	// Calculate exponential weighting for these three monitors
	measurement weightSP = expWeightMonVect(spCts);
	measurement weightOld = expWeightMonVect(oldCts);
	measurement weightBare = expWeightMonVect(bareCts);

	// Figure out when the dagger moves
	std::vector<double> dagSteps = dagDips(mcs1,mcs2);	
	
	// If there are no dagger movements, skip the run (it was probably killed due to bad fill) 
	if (dagSteps.end() != dagSteps.begin()) { 
		
		// Begin finding background counts
			// Need to modify when background actually happens:
			// For our past run cycle, this has been the last 50s of the run.
			// At some point we should also generate position-dependent backgrounds.
		//printf("Dagger moves at:\n!");
		//for (int j = 0; j < dagSteps.size(); j++) {
		//	printf("%f\n",*(dagSteps.end()-2));
		//}
		//double bkgStart;
		
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
		
		std::vector<input_t> bkgDagCts = mcs1->getCounts(
			[](input_t x)->input_t{return x;},
			[bkgStart,bkgEnd](input_t x)->bool{return ((x.ch == 1 || x.ch == 2) && x.realtime > bkgStart && x.realtime < bkgEnd);}
		);
		//std::vector<input_t> bkgDagCts = mcs1->getCoincCounts(
		//	[](input_t x)->input_t{return x;},
		//	[bkgStart, bkgEnd](input_t x)->bool{return (x.realtime > bkgStart && x.realtime < bkgEnd);}
		//);
		
		// Output a .csv file that we can use
		FILE* outfile;
		outfile = fopen("normNByDip.csv", "a");
		printf("Outputting data to normNByDip.csv\n!");
		
		// Loop through each dagger step and find when photon events hit
		//for(auto stepIt = dagSteps.begin()+1; stepIt < dagSteps.end(); stepIt++) {
		for(auto stepIt = dagSteps.begin()+1; stepIt < dagSteps.end(); stepIt++) {
			double startTime = *(stepIt-1);
			double endTime = *(stepIt);
			//fprintf(stderr, "\nFinding Start and Stop times for last dip for Field Scan!!!!!!!!!!!!!!!!!!!!\n");
			
			std::vector<input_t> dagCts = mcs1->getCounts(
				[](input_t x)->input_t{return x;},
				[startTime, endTime](input_t x)->bool{return ((x.ch == 1 || x.ch == 2) && x.realtime > startTime && x.realtime < endTime);}
			);
			
			if(stepIt == dagSteps.end()-1) { 
				endTime = bkgStart;
			}
			//if(stepIt == dagSteps.end()-1) {
			//	endTime = dagCts.back().realtime;
			//}
			
			// Added in retrieving coincidences to calculate step mean - more stable due to Bkg.
			std::vector<input_t> coinc = mcs1->getCoincCounts(
				[](input_t x)->input_t{return x;},
				[startTime, endTime](input_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
			);
			
			double stepMean = coinc.size() > 0 ?
				std::accumulate(coinc.begin(), coinc.end(), 0.0, [](double m, input_t x)->double{return m + x.realtime;}) / (double)coinc.size()
				: 0.0;
			
			// Deadtime Correction 
			double time = floor(dagCts.front().realtime);
			double counts = 0.0;
			double corr = 0.0;
			double deadTimeCounts = getDeadTimeCounts(mcs1, startTime, endTime);
			for(auto cIt = dagCts.begin(); cIt < dagCts.end(); cIt++) {
				if(floor(cIt->realtime) != time) {
					corr += (counts/(1.0-counts*10*NANOSECOND) - counts); // Deadtime correction
					counts = 0.0;
					time = floor(cIt->realtime);
				}
				counts += 1.0;
			}
	
			//printf("Data -  %d,%f,%f,%f,%f,%lu,%f,%lu,%f,%lu,%f,%f,%f,%f,%f,%f,%f,%f\n",
			printf("\nData -----  Run: %d\ntdTime: %f\nstartTime: %f\nendTime: %f\nstepMean: %f\ncoinc: %lu\ndeadTime: %f\ndagCts: %lu\ncorr: %f\nbkgDag: %lu\nbkgTime: %f\nSP: %f (%f)\nOLD: %f(%f)\nBARE: %f (%f)\nfillEnd: %f\n",
				   mcs1->getRunNo(), tdTime, startTime, endTime, stepMean, coinc.size(), deadTimeCounts,
				   dagCts.size(), corr, bkgDagCts.size(), bkgEnd - bkgStart,
				   weightSP.val, weightSP.err, weightOld.val, weightOld.err,
				   weightBare.val, weightBare.err, fillEnd);
			fprintf(outfile,"%d,%f,%f,%f,%f,%lu,%f,%lu,%f,%lu,%f,%f,%f,%f,%f,%f,%f,%f\n",
				   mcs1->getRunNo(), tdTime, startTime, endTime, stepMean, coinc.size(), deadTimeCounts,
				   dagCts.size(), corr, bkgDagCts.size(), bkgEnd - bkgStart,
				   weightSP.val, weightSP.err, weightOld.val, weightOld.err,
				   weightBare.val, weightBare.err, fillEnd);
		}
		fclose(outfile);
		fitFill(mcs1, mcs2);
	}
}*/

/*void longRunBkg(Run* run) {
 * // SHITTY VERSION CONTROL HERE
 * 
	std::vector<double> dagSteps = dagDips(run, NULL);
	if(dagSteps.empty()) {
		return;
	}
	double holdEnd = dagSteps.front()-10;
	double holdStart = dagSteps.front()-1010;
	double firstPeakStart = dagSteps.size() > 2 ? dagSteps[0] : 0.0;
	double firstPeakEnd = dagSteps.size() > 2 ? dagSteps[1] : 0.0;
	double endBkgStart = *(dagSteps.end()-2) + 100;
	double endBkgEnd = *(dagSteps.end()-1);
	
	unsigned long nHold[3];
	unsigned long nFirst[3];
	unsigned long nEnd[3];
	
	std::vector<input_t> ctsD1;
	std::vector<input_t> ctsD2;
	std::vector<input_t> ctsCo;

	std::vector<input_t> bkgDagCts = mcs1->getCoincCounts(
		[](input_t x)->input_t{return x;},
		[bkgStart, bkgEnd](input_t x)->bool{return (x.realtime > bkgStart && x.realtime < bkgEnd);}
	);
	// Find the end/start of the hold, taking out to 1200.1 to count 1199
	ctsD1 = run->getCounts(
		[](input_t x)->input_t{return x;},
		[holdEnd, holdStart](input_t x)->bool{return x.ch == 1 && x.realtime > holdStart && x.realtime < holdEnd;}
	);
	ctsD2 = run->getCounts(
		[](input_t x)->input_t{return x;},
		[holdEnd, holdStart](input_t x)->bool{return x.ch == 2 && x.realtime > holdStart && x.realtime < holdEnd;}
	);
	ctsCo = run->getCoincCounts(
		[](input_t x)->input_t{return x;},
		[holdEnd, holdStart](input_t x)->bool{return x.realtime > holdStart && x.realtime < holdEnd;}
	);
	nHold[0] = ctsD1.size();
	nHold[1] = ctsD2.size();
	nHold[2] = ctsCo.size();
		
	// find the end/start of the first peak, again taking out to 1200.1
	ctsD1 = run->getCounts(
		[](input_t x)->input_t{return x;},
		[firstPeakEnd, firstPeakStart](input_t x)->bool{return x.ch == 1 && x.realtime > firstPeakStart && x.realtime < firstPeakEnd;}
	);
	ctsD2 = run->getCounts(
		[](input_t x)->input_t{return x;},
		[firstPeakEnd, firstPeakStart](input_t x)->bool{return x.ch == 2 && x.realtime > firstPeakStart && x.realtime < firstPeakEnd;}
	);
	ctsCo = run->getCoincCounts(
		[](input_t x)->input_t{return x;},
		[firstPeakEnd, firstPeakStart](input_t x)->bool{return x.realtime > firstPeakStart && x.realtime < firstPeakEnd;}
	);
	nFirst[0] = ctsD1.size();
	nFirst[1] = ctsD2.size();
	nFirst[2] = ctsCo.size();
	
	// find the end of the background run 
	ctsD1 = run->getCounts(
		[](input_t x)->input_t{return x;},
		[endBkgEnd, endBkgStart](input_t x)->bool{return x.ch == 1 && x.realtime > endBkgStart && x.realtime < endBkgEnd;}
	); 
	
	ctsD2 = run->getCounts(
		[](input_t x)->input_t{return x;},
		[endBkgEnd, endBkgStart](input_t x)->bool{return x.ch == 2 && x.realtime > endBkgStart && x.realtime < endBkgEnd;}
	);
	ctsCo = run->getCoincCounts(
		[](input_t x)->input_t{return x;},
		[endBkgEnd, endBkgStart](input_t x)->bool{return x.realtime > endBkgStart && x.realtime < endBkgEnd;}
	);
	nEnd[0] = ctsD1.size();
	nEnd[1] = ctsD2.size();
	nEnd[2] = ctsCo.size();
	
	double lastCount = 0.0;
	if(ctsD1.size() > 0 && ctsD2.size() > 0) {
		lastCount = ctsD1.back().realtime > ctsD2.back().realtime ? ctsD1.back().realtime : ctsD2.back().realtime;
	}
	
	printf("Data - %d,%lu,%lu,%lu,%f,%lu,%lu,%lu,%f,%lu,%lu,%lu,%f\n", run->getRunNo(),
		  nHold[0], nHold[1], nHold[2], holdEnd-holdStart,
		  nFirst[0], nFirst[1], nFirst[2], firstPeakEnd-firstPeakStart,
		  nEnd[0], nEnd[1], nEnd[2], lastCount-endBkgStart);
	
	return;
}
*/

/* Object for fitting difference of fitted histograms */
/*void getDiffFitHist(TH1* hist, TF1* fit, double& neg, double& pos) {
	
	pos = 0.0;
	neg = 0.0;
	double diff = 0.0;
	
	for(int i = 1; i < hist->getNBinsX() +1; i++) {
		diff = hist->GetBinContent(i) - fit->Eval(hist->GetBinCenter(i));
		pos = diff > pos ? diff : pos;
		neg = diff < neg ? diff : neg;
	}
}*/

/* Object for fitting average offset of fitted histograms */
/*void getDiffAvgHist(TH1* hist, double avg, double& neg, double& pos) {
	
	pos = 0.0;
	neg = 0.0;
	double diff = 0.0;
	
	for(int i = 1; i < hist->getNBinsX() +1; i++) {
		diff = hist->GetBinContent(i) - avg;
		pos = diff > pos ? diff : pos;
		neg = diff < neg ? diff : neg;
	}
}*/
