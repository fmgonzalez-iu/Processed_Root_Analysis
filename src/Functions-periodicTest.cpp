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

	These functions are something that have to do with periodic stuff.
	I need more time to understand these.
	* 
------------------------------------------------------------------------------*/

/* Create two upsilon (frequency space?) histograms from our runs */
void protheroePeriodicDagStep(Run* runMCS1, Run* runMCS2, TH1D* histA, TH1D* histB) {
	// initialize variables
	int numCts = 100.0;
	double upsilonA = 0.0;
	double upsilonB = 0.0;
	double freq = 20003.0205;
	double stepStart = runMCS1->getTagBitEvt(1<<9, 160-2, 0);
	double stepEnd = runMCS1->getTagBitEvt(1<<9, stepStart + 0.25, 1);
	double start = stepEnd + 1.0;

	std::vector<input_t> ctsA = runMCS2->getCounts(
		[](input_t x)->input_t{return x;},
		[stepStart, stepEnd](input_t x)->bool{return x.ch == 1 && x.realtime > stepStart && x.realtime < stepEnd;}
	);
	std::vector<input_t> ctsB = runMCS2->getCounts(
		[](input_t x)->input_t{return x;},
		[stepStart, stepEnd](input_t x)->bool{return x.ch == 1 && x.realtime > stepStart && x.realtime < stepEnd;}
	);
	
	printf("%d %ld %ld\n", runMCS1->getRunNo(), ctsA.size(), ctsB.size());
	
	if(ctsA.size() < numCts || ctsB.size() < numCts) {
		printf("Not enough counts!\n");
		return;
	}

	// Changing our two vectors from counts to frequency space (?)
	printf("transforming\n");
	std::vector<double> phiA;
	std::transform(ctsA.end()-numCts, ctsA.end(), back_inserter(phiA), [freq](input_t x)->double{return fmod(x.realtime, 1.0/freq)/(1.0/freq);});
	std::vector<double> phiB;
	std::transform(ctsB.end()-numCts, ctsB.end(), back_inserter(phiB), [freq](input_t x)->double{return fmod(x.realtime, 1.0/freq)/(1.0/freq);});
	
	double numA = (double)phiA.size();
	double numB = (double)phiB.size();
	
	// summing our two vectors (again in possibly frequency space?)
	printf("summing\n");
	for(auto ii = phiA.begin(); ii < phiA.end()-1; ii++) {
		for(auto jj = ii+1; jj < phiA.end(); jj++) {
			upsilonA += 2.0/((0.5-fabs(fabs(*ii-*jj)-0.5)+1.0/numA)*(numA*(numA-1.0)));
		}
	}
	for(auto ii = phiB.begin(); ii < phiB.end()-1; ii++) {
		for(auto jj = ii+1; jj < phiB.end(); jj++) {
			upsilonB += 2.0/((0.5-fabs(fabs(*ii-*jj)-0.5)+1.0/numB)*(numB*(numB-1.0)));
		}
	}
	printf("Finished Sum\n");
	histA->Fill(upsilonA);
	histB->Fill(upsilonB);
	
	printf("Data - %d,%f,%f,%ld,%ld\n", runMCS1->getRunNo(), upsilonA, upsilonB, phiA.size(), phiB.size());
}


/* creates a single (frequency space) transform or something? 
 * fun fact: a search for Protheroe comes up with a British song from 1974 called "Pinball". */
void protheroePeriodicTest(Run* run, TH1D* hist) {
	// declare variables 
	int i = 0;
	int maxIA = 0;
	int maxIB = 0;
	double freq = 20003.0205;
	double maxA = 0.0;
	double maxB = 0.0;
	double minDelta = 10.0;
	std::vector<double> upsilonAs;
	std::vector<double> upsilonBs;
	std::vector<double> coincTimes;
	
	// create vector of counts. Take counts out to 1200.1 so 1199 gets counted
	std::vector<input_t> ctsA = run->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 1 && x.realtime > 200.0 && x.realtime < 1200.1;}
	);
	std::vector<input_t> ctsB = run->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 2 && x.realtime > 200.0 && x.realtime < 1200.1;}
	);
	std::vector<coinc_t> coinc = run->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[](coinc_t x)->bool{return x.realtime > 200.0 && x.realtime < 1200.0;}
	);
	
	std::transform(coinc.begin(), coinc.end(), back_inserter(coincTimes), [](coinc_t x)->double{return floor(x.realtime);});
	printf("%ld\n", coincTimes.size());
	
	auto prevIt = ctsA.begin();
	for(auto ii = ctsA.begin(); ii < ctsA.end(); ii++) {
		if(floor(prevIt->realtime) != floor(ii->realtime)) {
			// Check to see if we're in a coincidence period
			if(std::find(coincTimes.begin(), coincTimes.end(), floor(prevIt->realtime)) == coincTimes.end()) {
				std::vector<double> phiA;
				std::transform(prevIt, ii, back_inserter(phiA), [freq](input_t x)->double{return fmod(x.realtime, 1.0/freq)/(1.0/freq);});
				double upsilonA = 0.0;
				if(phiA.size() >= 300) {
					auto end = phiA.begin() + 300;
					int numA = std::distance(phiA.begin(), end);
					// I renamed the iterators, because ii and jj previously had the same name?
					for(auto jj = phiA.begin(); jj < end - 1; jj++) {
						for(auto kk = jj+1; kk < end; kk++) {
							upsilonA += 2.0/((0.5-fabs(fabs(*jj-*kk)-0.5)+1.0/numA)*(numA*(numA-1.0)));
						}
					}
					upsilonAs.push_back(upsilonA);
					hist->Fill(upsilonA);
					maxIA = upsilonA > maxA ? floor(prevIt->realtime) : maxIA;
					maxA = upsilonA > maxA ? upsilonA : maxA;
				}
				prevIt = ii;
			}
			else {
				prevIt = ii;
			}
		}
	}
	double avgA = std::accumulate(upsilonAs.begin(), upsilonAs.end(), 0.0) / upsilonAs.size();
	double stdDevA = sqrt(std::accumulate(upsilonAs.begin(), upsilonAs.end(), 0.0, [avgA, &upsilonAs](double s, double u)->double{return s + (1.0/upsilonAs.size())*(u-avgA)*(u-avgA);}));
	printf("%d,%.18f,%.18f,%.18f,%d,%.18f,%.18f\n", run->getRunNo(), avgA, maxA, stdDevA, maxIA, (maxA-avgA)/stdDevA, minDelta*(1.0/freq)/NANOSECOND);
}

/* Find the rayleigh background */
void rayleighPeriodicTestBkg(Run* run) {
	double freq = 20003.02232;
	std::vector<input_t> ctsA = run->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 1 && x.realtime > 200.0 && x.realtime < 1200.0;}
	);
	std::vector<input_t> ctsB = run->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 2 && x.realtime > 200.0 && x.realtime < 1200.0;}
	);
	std::vector<coinc_t> coinc = run->getCoincCounts(
		[](coinc_t x)->coinc_t{return x;},
		[](coinc_t x)->bool{return x.realtime > 200.0 && x.realtime < 1200.0;}
	);
	std::vector<double> coincTimes;
	std::transform(coinc.begin(), coinc.end(), back_inserter(coincTimes), [](coinc_t x)->double{return floor(x.realtime);});
	printf("%ld\n", coincTimes.size());
	std::vector<double> cosSums;
	std::vector<double> sinSums;
	for(int i = 0; i < 20; i++) {
		sinSums.push_back(0.0);
		cosSums.push_back(0.0);
	}
	int num = 0;
	double maxA = 0.0;
	double maxATime = 0.0;
	double maxBTime = 0.0;
	double maxB = 0.0;
	double avgA = 0.0;
	int numA = 0;
	double avgB = 0.0;
	int numB = 0;
	for(auto it = ctsA.begin()+1; it < ctsA.end(); it++) {
		if(floor(it->realtime) != floor((it-1)->realtime)) {
			double z2_20 = 2.0*(
				std::accumulate(sinSums.begin(), sinSums.end(), 0.0, [](double s, double x){return s+x*x;})
				+std::accumulate(cosSums.begin(), cosSums.end(), 0.0, [](double s, double x){return s+x*x;})
			) / float(num);
			if(std::find(coincTimes.begin(), coincTimes.end(), floor((it-1)->realtime)) == coincTimes.end()) {
				avgA += z2_20;
				maxATime = z2_20 > maxA ? floor((it-1)->realtime) : maxATime;
				maxA = z2_20 > maxA ? z2_20 : maxA;
				numA += 1;
			}
			num = 0;
			for(int i = 0; i < 20; i++) {
				sinSums[i] = 0.0;
				cosSums[i] = 0.0;
			}
		}
		num += 1;
		for(int i = 1; i <= 20; i++) {
			sinSums[i-1] += sin(i*2.0*M_PI*freq*it->realtime);
			cosSums[i-1] += cos(i*2.0*M_PI*freq*it->realtime);
		}
	}
	for(int i = 0; i < 20; i++) {
		sinSums[i] = 0.0;
		cosSums[i] = 0.0;
	}
	num = 0;
	for(auto it = ctsB.begin()+1; it < ctsB.end(); it++) {
		if(floor(it->realtime) != floor((it-1)->realtime)) {
			double z2_20 = 2.0*(
				std::accumulate(sinSums.begin(), sinSums.end(), 0.0, [](double s, double x){return s+x*x;})
				+std::accumulate(cosSums.begin(), cosSums.end(), 0.0, [](double s, double x){return s+x*x;})
			) / float(num);
			if(std::find(coincTimes.begin(), coincTimes.end(), floor((it-1)->realtime)) == coincTimes.end()) {
				avgB += z2_20;
				maxBTime = z2_20 > maxB ? floor((it-1)->realtime) : maxBTime;
				maxB = z2_20 > maxB ? z2_20 : maxB;
				numB += 1;
			}
			num = 0;
			for(int i = 0; i < 20; i++) {
				sinSums[i] = 0.0;
				cosSums[i] = 0.0;
			}
		}
		num += 1;
		for(int i = 1; i <= 20; i++) {
			sinSums[i-1] += sin(i*2.0*M_PI*freq*it->realtime);
			cosSums[i-1] += cos(i*2.0*M_PI*freq*it->realtime);
		}
	}
	printf("Data - %d,%f,%f,%f,%f,%f,%f\n", run->getRunNo(), avgA / numA, maxA, maxATime, avgB / numB, maxB, maxBTime);
}

/* test of the raleigh function */ 
void rayleighPeriodicTest(Run* run) {
	double freq = 20003.02232;
	double z2_20_A = 0.0;
	double z2_20_B = 0.0;
	std::vector<input_t> ctsA = run->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 1 && x.realtime > 200.0 && x.realtime < 1200.0;}
	);
	std::vector<input_t> ctsB = run->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 2 && x.realtime > 200.0 && x.realtime < 1200.0;}
	);
	
	for(int i = 1; i <= 20; i++) {
		z2_20_A += 2.0 * (
			pow(std::accumulate(ctsA.begin(), ctsA.end(), 0.0, [i, freq](double r, input_t x)->double{
			return r + sin(i*2.0*M_PI*freq*x.realtime);}), 2.0)
			+ pow(std::accumulate(ctsA.begin(), ctsA.end(), 0.0, [i, freq](double r, input_t x)->double{
			return r + cos(i*2.0*M_PI*freq*x.realtime);}), 2.0)
			) / ctsA.size();
	}
	for(int i = 1; i <= 20; i++) {
		z2_20_B += 2.0 * (
			pow(std::accumulate(ctsB.begin(), ctsB.end(), 0.0, [i, freq](double r, input_t x)->double{
			return r + sin(i*2.0*M_PI*freq*x.realtime);}), 2.0)
			+ pow(std::accumulate(ctsB.begin(), ctsB.end(), 0.0, [i, freq](double r, input_t x)->double{
			return r + cos(i*2.0*M_PI*freq*x.realtime);}), 2.0)
			) / ctsB.size();
	}
	printf("Data - %d,%.17f,%.17f,%ld,%ld\n", run->getRunNo(), z2_20_A, z2_20_B, ctsA.size(), ctsB.size());
}
