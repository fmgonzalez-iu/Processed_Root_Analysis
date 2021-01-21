#include "../inc/Functions.hpp"
#include "../inc/ExpFill.hpp"
#include "TCanvas.h"
#include "TGraph.h"

/* define constants we need for later */
#define NANOSECOND .000000001
#define INVKAPPA 1.0/70.0 // inverse of trap draining time hardcode

/*----------------------------------------------------------------------------
	Author: Frank M. Gonzalez
	
	This file contains functions to quantify the exponential weighting.
	The number of UCN trapped in the experiment can be found by:
	 
	N_{trapped} = N_{max} ( 1 - Exp[-t/\kappa])
	
	* \kappa is  the saturation time constant in the trap.
	
	For a constant flux \Phi, we can integrate to get the density:
	
	N_{trapped} = \gamma \Int[\Phi Exp[(t-T_{end})/\kappa]]
	
	* \gamma is a proportionality constant
	* \Int[\Phi] is replaced by \Sum[(1/\kappa)] over times t_i of individual pulses.
	* 
	N_{trapped} = \gamma \Phi \Sum[(1/\kappa) Exp[(t_i-T_{end}) /\kappa]
	
	From these we can get a weighted monitor signal.
	The saturation time of the trap, \kappa, is measured as 70.0s, according to 2014 fits.
------------------------------------------------------------------------------*/

/* Set the exponential weighting.
 * Input is a histogram of monitor counts and the end of loading */
measurement expWeightMon(TH1D* mon, double end) {
	
	// measurement consists of two things -- value and error.
	measurement weight = {0.0, 0.0};
	
	// loop across all bins and fill in our weighting function, using the sum defined above
	int i = 0;
	int nBinsX = mon->GetNbinsX();
	for(i = 1; i <= nBinsX; i++) {
		weight.val += mon->GetBinContent(i)*INVKAPPA*exp(INVKAPPA*(mon->GetBinCenter(i)-end));
		weight.err += mon->GetBinContent(i)*INVKAPPA*INVKAPPA*exp(2*INVKAPPA*(mon->GetBinCenter(i)-end));	
	}
	weight.err = sqrt(weight.err);
	return weight;
}

/* Another weight  on this monitor monitor, but this time loading a vector of our counts */
measurement expWeightMonVect(std::vector<input_t> &cts) {
	
	// measurement consists of two things -- value and error.
	measurement weight = {0.0, 0.0};

	// We start at a rate of 0 at the front of our vector of counts
	int rate = 0;
	double prevTime;
	if (cts.size() != 0) {
		prevTime = floor(cts.front().realtime);
	} else { return weight; }
	double dtCts = 0.0;
	// loop across all bins and fill in our weighting function, using the sum defined above
	for(auto it = cts.begin(); it < cts.end(); it++) {
		weight.val += INVKAPPA*exp(INVKAPPA * it->realtime);
		weight.err += INVKAPPA*INVKAPPA*exp(2.0 * INVKAPPA * it->realtime);
		
		// Also account for deadtime in the detector. Take the inv. kappa of the 
		if(floor(it->realtime) > prevTime || it == cts.end()-1) {
			dtCts = (rate/(1.0-rate*20*NANOSECOND) - rate); // arbitrary 20ns deadtime
			weight.val += dtCts * INVKAPPA * exp(INVKAPPA * (prevTime + 0.5));
			weight.err += dtCts * INVKAPPA*INVKAPPA * exp(2.0 * INVKAPPA * (prevTime+0.5));
			
			// Set prevTime to the next value
			prevTime = floor(it->realtime);
			rate = 0;
		}
		rate +=1;
	}
	weight.err = sqrt(weight.err);
	return weight;
}

/* Another weight  on this monitor monitor, but this time loading a vector of our counts */
measurement expWeightMonVect(std::vector<input_t> &cts, double kappa) {
	
	// measurement consists of two things -- value and error.
	measurement weight = {0.0, 0.0};
	
	// First make sure we're given a valid time constant for this run
	double timeC = 1.0/kappa;
	if (timeC < 0.0) {
		return weight;
	}
	
	// We start at a rate of 0 at the front of our vector of counts
	int rate = 0;
	double prevTime;
	if (cts.size() != 0) {
		prevTime = floor(cts.front().realtime);
	} else { return weight; }
	
	double dtCts = 0.0;
	// loop across all bins and fill in our weighting function, using the sum defined above
	for(auto it = cts.begin(); it < cts.end(); it++) {
		weight.val += timeC*exp(timeC * it->realtime);
		weight.err += timeC*timeC*exp(2.0 * timeC * it->realtime);
		
		// If our counts vector bin occured after the previous time, reweight.
		// We're doing an average across whatever our previous time was.
		if(floor(it->realtime) > prevTime || it == cts.end()-1) {
			dtCts = (rate/(1.0-rate*20*NANOSECOND) - rate); // arbitrary 20ns deadtime
			weight.val += dtCts * timeC * exp(timeC * prevTime + 0.5);
			weight.err += dtCts * timeC*timeC * exp(2.0 * timeC * (prevTime+0.5));
			
			// Set prevTime to the next value
			prevTime = floor(it->realtime);
			rate = 0;
		}
		rate +=1;
	}
	weight.err = sqrt(weight.err);
	return weight;
}

/* Another weight but now incorporating variable deadtime */
measurement expWeightMonVectDT(std::vector<input_t> &cts, double kappa, double deadtime, double end) {
	
	// measurement consists of two things -- value and error.
	measurement weight = {0.0, 0.0};
	
	// We start at a rate of 0 at the front of our vector of counts
	double rate = 0.0;
	//double prevTime;
	double dtS = 0.1; // want deadtime averaged every 0.1s
	int time;
	if (cts.size() != 0) { // This is how we check dead time
		time = floor(cts.front().realtime);
	} else { return weight; }
	
	
	//int time=floor(cts.front().realtime / dtS);
	// First make sure we're given a valid time constant for this run
	if ((-1e-7 < kappa) && (kappa < 1e-7)) { // If kappa is zero, just sum raw cts 
		// But we're still doing a deadtime correction	
		double rate = 0.0;
		for (auto it = cts.begin(); it < cts.end(); it++) {
			weight.val += 1.0; // Everything is weighted as 1 though 
			weight.err += 1.0;
			if(floor(it->realtime / dtS) != time) {
				double corr = (rate/(1.0-rate*deadtime / dtS) - rate); 
				weight.val += corr; // So is the deadtime
				weight.err += corr;
				rate = 0.0;
				time = floor(it->realtime / dtS);
			}
			rate += 1.0;
		}
		weight.err = sqrt(weight.err);
		return weight; // We can end here if there's no time constant
	}
	
	double timeC = 1.0/kappa; // If we have a fitted constant, go here
	if (timeC < 0.0) {
		return weight;
	}
	
	//double dtCts = 0.0;
	//double rate = 0.0;
	// loop across all bins and fill in our weighting function, using the sum defined above
	for(auto it = cts.begin(); it < cts.end(); it++) {
		weight.val += timeC*exp(timeC * (it->realtime - end)); // Note that in the other functions, end is not defined
		weight.err += timeC*timeC*exp(2.0 * timeC * (it->realtime - end)); // This just acts as a constant so it doesn't actually affect anything
		
		// If our counts vector bin occured after the previous time, reweight.
		// We're doing an average across whatever our previous time was.
		if(floor(it->realtime / dtS) > time || it == cts.end()-1) {
			double corr = (rate/(1.0-rate*deadtime / dtS) - rate); // correct rate for deadtime
			weight.val += corr * timeC * exp(timeC * ((time + 0.5)*dtS - end));
			weight.err += corr * timeC*timeC * exp(2.0 * timeC * ((time+0.5)*dtS - end));
			
			// Set time to the next value
			time = floor(it->realtime / dtS);
			rate = 0;
		}
		rate +=1;
	}
	weight.err = sqrt(weight.err);
	return weight;
}

measurement expWeightMonVectBkg(std::vector<input_t> &cts, double kappa, double deadtime, double bkg, double end) {
	
	// measurement consists of two things -- value and error.
	measurement weight = {0.0, 0.0};
	
	// We start at a rate of 0 at the front of our vector of counts
	double rate = 0.0;
	//double prevTime;
	double dtS = 0.1; // want deadtime averaged every 0.1s
	int time;
	if (cts.size() != 0) { // This is how we check dead time
		time = floor(cts.front().realtime);
	} else { return weight; }
	
	
	//int time=floor(cts.front().realtime / dtS);
	// First make sure we're given a valid time constant for this run
	if ((-1e-7 < kappa) && (kappa < 1e-7)) { // If kappa is zero, just sum raw cts 
		// But we're still doing a deadtime correction	
		double rate = 0.0;
		for (auto it = cts.begin(); it < cts.end(); it++) {
			weight.val += 1.0; // Everything is weighted as 1 though 
			weight.err += 1.0;
			if(floor(it->realtime / dtS) != time) {
				double corr = (rate/(1.0-rate*deadtime / dtS) - rate); 
				weight.val += corr; // So is the deadtime
				weight.err += corr;
				weight.val -= bkg*dtS;  // Subtract off background
				weight.err += bkg*dtS;
				rate = 0.0;
				time = floor(it->realtime / dtS);
			}
			rate += 1.0;
		}
		weight.err = sqrt(weight.err);
		return weight; // We can end here if there's no time constant
	}
	
	double timeC = 1.0/kappa; // If we have a fitted constant, go here
	if (timeC < 0.0) {
		return weight;
	}
	
	//double dtCts = 0.0;
	//double rate = 0.0;
	// loop across all bins and fill in our weighting function, using the sum defined above
	for(auto it = cts.begin(); it < cts.end(); it++) {
		weight.val += timeC*exp(timeC * (it->realtime - end)); // Note that in the other functions, end is not defined
		weight.err += timeC*timeC*exp(2.0 * timeC * (it->realtime - end)); // This just acts as a constant so it doesn't actually affect anything
		
		// If our counts vector bin occured after the previous time, reweight.
		// We're doing an average across whatever our previous time was.
		if(floor(it->realtime / dtS) > time || it == cts.end()-1) {
			double corr = (rate/(1.0-rate*deadtime / dtS) - rate); // correct rate for deadtime
			weight.val += corr * timeC * exp(timeC * ((time + 0.5)*dtS - end));
			weight.err += corr * timeC*timeC * exp(2.0 * timeC * ((time+0.5)*dtS - end));
			
			double bCts = bkg * dtS;
			weight.val -= bCts * timeC * exp(timeC * ((time - 0.5)*dtS - end));
			weight.err += bCts * timeC*timeC * exp(2.0 * timeC * ((time+0.5)*dtS - end));
			
			// Set time to the next value
			time = floor(it->realtime / dtS);
			rate = 0;
		}
		rate +=1;
	}
	weight.err = sqrt(weight.err);
	return weight;
}

/*  Maybe try using fit phi values from expfill to get monitor values? */
measurement expWeightMonVectPhi(std::vector<input_t> &cts, double kappa) {
	
	// measurement consists of two things -- value and error.
	measurement weight = {0.0, 0.0};
	double timeC = 1.0/kappa;
/*double ExpFill::operator() (double *x, double *p) {
	double phi = 0.0;
	int i = 0;
	for(i = 0; i < hgxHits.size(); i++) {
		float T_i = (hgxHits[i]+1.26); // hardcoded in time offset.
		if(x[0] < T_i) {
			continue;
		}
		// time constants RH mon for 2018.
		// use ExpFillFree to find new time constants.
		phi += p[i]
			*(1.0 - exp(-(x[0]-T_i)/2.20))
			*exp(-(x[0]-T_i)/21.38);
	}
	return phi;
}*/

	// We start at a rate of 0 at the front of our vector of counts
	int rate = 0;
	double prevTime = floor(cts.front().realtime);
	double dtCts = 0.0;
	// loop across all bins and fill in our weighting function, using the sum defined above
	for(auto it = cts.begin(); it < cts.end(); it++) {
		weight.val += timeC*exp(timeC * it->realtime);
		weight.err += timeC*timeC*exp(2.0 * timeC * it->realtime);
		
		// If our counts vector bin occured after the previous time, reweight.
		// We're doing an average across whatever our previous time was.
		if(floor(it->realtime) > prevTime || it == cts.end()-1) {
			dtCts = (rate/(1.0-rate*20*NANOSECOND) - rate); // arbitrary 20ns deadtime
			weight.val += dtCts * timeC * exp(timeC * prevTime + 0.5);
			weight.err += dtCts * timeC*timeC * exp(2.0 * timeC * (prevTime+0.5));
			
			// Set prevTime to the next value
			prevTime = floor(it->realtime);
			rate = 0;
		}
		rate +=1;
	}
	weight.err = sqrt(weight.err);
	return weight;
}
