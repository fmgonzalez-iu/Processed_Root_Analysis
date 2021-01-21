#include "../inc/ExpFill.hpp"

/*	------------------------------------------------------------------------------------------------
	Author: Frank M. Gonzalez
	
	This class characterizes the fill quality of the trap. 
	Arrival time is fit to a function:
	
	\Phi(t) = \Sum[\phi_i (1 - Exp[-(t-T_i)/\kappa_1]) * (Exp[-(t-T_i)/\kappa_2])]
	
	* \phi_i   = the strength of an individual pulse
	* T_i      = time offset of an individual pulse
	* \kappa_1 = saturation time constant
	* \kappa_2 = decay time constant
		
	For 2018 run cycle, time constants are chosen from a weekend's worth running (October 7):
	GV monitor (Detector 3)
	* OFFSET   = 1.47  +- (0.001)
	* \kappa_1 = 0.82  +- (0.007)
	* \kappa_2 = 20.47 +- (0.12)
	RH monitor (Detector 4)
	* OFFSET   = 1.26  +- (0.002)
	* \kappa_1 = 2.20  +- (0.009)
	* \kappa_2 = 21.38 +- (0.02)
	SP monitor (Detector 5)
	* OFFSET   = 2.73  +- (0.006)
	* \kappa_1 = 5.06  +- (0.46)
	* \kappa_2 = 11.44 +- (0.30)
	
	Note that the SP monitor has few enough counts that we can't really pin it down.
	Odd note that the RH monitor fits to a shorter "offset" than the GV.
	This is probably a fitting thing, with an emphasis on saturation of the monitor.
	------------------------------------------------------------------------------------------------	*/

/* initialize ExpFill() */
ExpFill::ExpFill() {
	numOffsets = 0;
	offsetT = 1.47; // GV monitor 2018
	kappa1  = 0.82;
	kappa2  = 20.47; 
}

//----------------------------------------------------------------------
/* set the time constants used for fitting */
void ExpFill::setOffsetT(double off) {
	offsetT = off;
}
double ExpFill::getOffsetT() {
	return offsetT;
}
/* filling time */
void ExpFill::setKappa1(double k1) {
	kappa1 = k1;
}
double ExpFill::getKappa1() {
	return kappa1;
}
/* draining time */
void ExpFill::setKappa2(double k2) {
	kappa2 = k2;
}
double ExpFill::getKappa2() {
	return kappa2;
}

//----------------------------------------------------------------------
/* add an offset (H-GX pulse) to ExpFill */
void ExpFill::addOffset(double off) {
	numOffsets += 1;
	hgxHits.push_back(off);
}
/* get the number of offsets (H-GX pulses) added to ExpFill */
int ExpFill::getNumOffsets() {
	return numOffsets;
}
//----------------------------------------------------------------------

/* actually fill the trap */
double ExpFill::operator() (double *x, double *p) {
	double phi = 0.0;
	int i = 0;
	for(i = 0; i < hgxHits.size(); i++) {
		float T_i = (hgxHits[i]+offsetT); // time offset.
		if(x[0] < T_i) {
			continue;
		}
	
		phi += p[i]
			*(1.0 - exp(-(x[0]-T_i)/kappa1))
			*exp(-(x[0]-T_i)/kappa2);
	}
	return phi;
}

ExpFill::~ExpFill() {

}
