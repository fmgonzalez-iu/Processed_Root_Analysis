#include "Run.hpp"
#include <numeric>

#pragma once

class ExpFill {
	public:
		ExpFill();
		~ExpFill();
		
		// Set/get H-GX pulses
		void addOffset(double off);
		int getNumOffsets();
		
		// Set time constants
		void setOffsetT(double ofT);
		void setKappa1(double k1);
		void setKappa2(double k2);
		
		// Get time constants
		double getOffsetT();
		double getKappa1();
		double getKappa2();
		
		//Old time constants -
		//Sat - 5.75 Decay - 10.66 Offset - 3.0 [ 2015 ? ]
		//Sat - 0.4  Decay - 24.75 Offset - 3.0 [ 2016-2017]
		double operator() (double *x, double *p);
	
	private:
		std::vector<double> hgxHits;
		int numOffsets;
		double offsetT = 1.47;
		double kappa1  = 0.82;
		double kappa2  = 20.47;
};

class ExpFillFree{
	public:
		ExpFillFree();
		~ExpFillFree();
		void addOffset(double off);
		
		int getNumOffsets();
		double getOffsetI(int i);
		double operator() (double* x, double* p);
	
	private:
		std::vector<double> hgxHits;
		int numOffsets;
};

/*class ExpFillFree {
	public:
		ExpFillFree(int num);
	
		//Old time constants -
		//Sat - 5.75 Decay - 10.66 Offset - 3.0
		double operator() (double *x, double *p);
	
	private:
		int nPulse;
};*/
