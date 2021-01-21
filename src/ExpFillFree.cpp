#include "../inc/ExpFill.hpp"

/*----------------------------------------------------------------------------
	Author: Frank M. Gonzalez
	
	This class characterizes the fill quality of the trap. 
	Arrival time is fit to a function:
	
	\Phi(t) = \Sum[\phi_i (1 - Exp[-(t-T_i)/\kappa_1]) * (Exp[-(t-T_i)/\kappa_2])]
	
	* \phi_i   = the strength of an individual pulse
	* T_i      = time offset of an individual pulse
	* \kappa_1 = saturation time constant
	* \kappa_2 = decay time constant
	
	In this class, we're assuming an ideal fill.
	Instead we're fitting the pulse offset and our time constants.
------------------------------------------------------------------------------*/

/* initialize ExpFillFree() */
ExpFillFree::ExpFillFree() {
	numOffsets = 0;
	if (hgxHits.size() != 0) {
		for (int i = 0; i < hgxHits.size(); i++) {
			hgxHits.pop_back();
		}
	}
}

/* add an offset (H-GX pulse) to ExpFill */
void ExpFillFree::addOffset(double off) {
	numOffsets += 1;
	hgxHits.push_back(off);
}
/* get the number of offsets (H-GX pulses) added to ExpFill */
int ExpFillFree::getNumOffsets() {
	return numOffsets;
}

double ExpFillFree::getOffsetI(int i) {
	return hgxHits[i];
}
/* actually fill the trap */
double ExpFillFree::operator() (double* x, double* param) {
	double phi = 0.0; // background
	double Ti = 0.0; // Initial time of hit
	//int toff = 0;
	//int npulses = hgxHits.size();
			
	// Now do an easy ~50 dimensional fit
	for(int j = 0; j < hgxHits.size(); j++) {
		// set time offset as an initial parameter
		Ti = hgxHits[j] + param[0];
		// x[0] represents the zeroth dimension of a ROOT histogram, in this case time. Other dims aren't used!
		// A word of warning --- ROOT fitting will look wonky if you have more than 48 fit elements.
		// I've checked this in Python, and the fitting is working as expected.
		if(x[0] > Ti) { 
			phi += param[j+3] * (1.0 - exp(-(x[0]-Ti)/param[1])) * exp(-(x[0]-Ti)/param[2]);
		}
		Ti = 0.0; // clear value of Ti
	}
	
	return phi;
}
ExpFillFree::~ExpFillFree() {

	//delete (*offset);
}
