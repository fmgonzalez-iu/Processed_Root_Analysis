#include "Run.hpp"
#include <numeric>

#pragma once

class DrainTime {
	public:
		DrainTime();
		~DrainTime();
		double operator() (double *x, double *p);
	
};

class MultiExpBkg {
	public: 
		MultiExpBkg(size_t len);
		~MultiExpBkg();
		double operator() (double *x, double *p);
		
		void setBackground(double r);
		void setTimeConsts(std::vector<double> tcGuess);
		void setScaleFactors(std::vector<double> tsGuess);
		void guessScales(double peak);
		
		size_t getNTimeConsts();
		double getBackground();
		std::vector<double> getTimeConsts;
						
	private:
		size_t nTC;
	
		std::vector<double> timeConsts;
		std::vector<double> tcScales;
		double background;
		
};

class FillTime {
	public:
		FillTime();
		~FillTime();
		double operator() (double *x, double *p);
		
};
