#include "Run.hpp"
#include <numeric> 

#pragma once
#define NANOSECOND .000000001

double getDeadTimeCoinc(std::vector<coinc_t> &cts, double start, double end);
double getDeadTimeSing(std::vector<input_t> &cts, double start, double end);

class Pileup {
	
	public:
		Pileup(Run* run, double bkg, double bkgC, double start, double end);
		~Pileup();
		
		
		void LoadCoincStatsVec(std::vector<coinc_t> coin_vec);
		void LoadSingStatsVec(std::vector<input_t> pmt1, std::vector<input_t> pmt2);
		void LoadSingStatsVec(std::vector<input_t> pmt);
		void LoadSingStatsVec(std::vector<coinc_t> pmt);
		
		
		void SetBackground(double bkgRate);
		double GetBackground();
		
		void SetCoincBackground(double bkgRate);
		double GetCoincBackground();
		
		void CalculateTau();
		double GetTau();
		
		double GetAmplitude();
		
		double CalculatePileup(std::vector<coinc_t> coinc);
		
	private:
		
		void LoadCoincStats(Run* run);
		void LoadSingStats(Run* run);
		
		double ph_f; // Foreground Counts
		double ph_b; // Background Counts in Window
		double coinc; // Coincidences in run
		double coincB; // Background coincidences in window
		
		double dt_mu; // Average coincidence length
		double N_mu ; // Average coincidence number of photons
		
		double N_max; // Maximum number of photons in coincidence events
		
		int peThresh; // Photon threshold
		
		double bkgRate;
		double start;
		double end;
		double tau;
		
		std::vector<double> logGamma;
	
};

