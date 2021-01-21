#include "../inc/DrainTime.hpp"

/*----------------------------------------------------------------------------
	Author: Frank M. Gonzalez
	
	These classes characterize time contants of a given run.
	We can feed this through ROOT's fitting algorithm.
	* 
	MultiExpBkg is a way to find backgrounds by fitting to a multi-exponential
	* Assume the background is constant + sum of exponentials
	* Overprogrammed this to pull things out the object that might not 
	* really be needed
	DrainTime is the draining time for a given dagger movement:
	* Linear for some time
	* Exponential decay afterwards
	* Normalized so probability is 1.
	FillTime is the filling time constant for a single exponential.
	* 
	* 
------------------------------------------------------------------------------*/

//----------------------------------------------------------------------
// MultiExp Background Calculator
//----------------------------------------------------------------------
/*initiator & destructor */
MultiExpBkg::MultiExpBkg(size_t len) { 
	// Initialize a MultiExpBkg object with an initial size
	nTC = len;
	// Clear the time constant and scaling factors (safety)
	if (timeConsts.size() != 0) {
		for(size_t ii; ii < timeConsts.size(); ii++) {
			timeConsts.pop_back();
		}
	}
	if (tcScales.size() != 0) {
		for(size_t ii; ii< tcScales.size(); ii++) {
			tcScales.pop_back();
		}
	}
	// Now initialize the object
	for (size_t ii; ii < len; ii++){
		// Initialize at TC = 10s, scaled to 1000 Hz
		timeConsts.push_back(10.);
		tcScales.push_back(1000.);
	}
	// Initialize background at zero
	background = 0.;
}

MultiExpBkg::~MultiExpBkg() { }

// Set return private objects
size_t MultiExpBkg::getNTimeConsts() {
	return nTC;
}

void MultiExpBkg::setBackground(double r) {
	// Background will be a rate
	background = r;
}

double MultiExpBkg::getBackground() {
	return background;
}

void MultiExpBkg::setTimeConsts(std::vector<double> tcGuess) {
	if (tcGuess.size() == nTC) {
		for (size_t ii; ii < nTC; ii++) {
			timeConsts.at(ii) = tcGuess.at(ii);
		}
	}		
}

void MultiExpBkg::guessScales(double peak) {
	// We're going to guess that the N time constants start at the same 
	// height. This height is going to be the peak rate/N.
	for (size_t ii; ii < nTC; ii++) {
		tcScales.at(ii) = peak / ((double)nTC);
	}
}

void MultiExpBkg::setScaleFactors(std::vector<double> tsGuess) {
	// If we know the time constant (from e.g. fitting). 
	if (tsGuess.size() == nTC) {
		for (size_t ii; ii < nTC; ii++) {
			tcScales.at(ii) = tsGuess.at(ii);
		}
	}
}
	
// For actual ROOT fitting, we need to have this object with an operator()
double MultiExpBkg::operator() (double* x, double* param) {
	// x[0] is the zeroth dimension of a ROOT histogram, i.e. time. 
	// Since this fits a 1D hist, we don't use the other dims.
	/*if (sizeof(param) != 2*nTC + 1) { 
		printf("ERROR! Unable to fit to misshapen parameter!\n");
		return -1;
	}*/
	double rate = param[0]; // First is the background
	for (int ii = 0; ii < nTC; ii++) {
		rate += param[2*ii+1] * exp(-x[0]/param[2*(ii+1)]);
	}
	return rate;
}

//----------------------------------------------------------------------
// Drain Time
//----------------------------------------------------------------------
/* initiator & destructor */
DrainTime::DrainTime() { }

DrainTime::~DrainTime() { }

/* DrainTime function takes two parameters, a linear cutoff and an
 * exponential time constant.
 * It returns a normalized probability distribution. */
double DrainTime::operator() (double* x, double* param) {
	
	double prob = 0.0; 
	
	// x[0] represents the zeroth dimension of a ROOT histogram, in this case time. Other dims aren't used!	
	if (x[0] < param[0]) { // Linear fill
		prob = x[0] /(param[0] * param[1] * (param[0] / (2*param[1] + 1)));
	} else { // Exponential part.
		prob = exp(-(x[0]-param[0])/param[1]) /(param[1] * (param[0] / (2*param[1] + 1)));
	}
	
	return prob;
	
}

//----------------------------------------------------------------------
// Fill Time
//----------------------------------------------------------------------
/* initiator & destructor */
FillTime::FillTime() { }

FillTime::~FillTime() { }

double FillTime::operator() (double* x, double* param) {

	// Simple decaying exponential with parameter to scale. 
	double Ntr = param[0] * (1.0 - exp(-x[0]/param[1]));
	
	return Ntr;
}
