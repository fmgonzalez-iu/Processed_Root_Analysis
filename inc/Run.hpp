#include <vector>
#include <algorithm>
#include <string>
#include <math.h>
#include <fstream>
#include <iterator>
#include <functional>
#include <numeric>
#include "stdio.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TString.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"

/* "#pragma once" tells the compiler to only compile included files once
 * to prevent multiple locations for Run.hpp appearing */
 
#pragma once
#define RATEWINDOW 0.1 // Tracking rate for singles/coincidences on this window (s)
#define NANOSECOND .000000001

/* Create a structure that contains the input from our runs: 
 * the tvc*/
struct input_t {
	// From Data
	unsigned long time;
	double realtime;
	int ch;
	int tag;
	// Other, maybe un-initialized stuff
	double deadtime; // deadtime in s
	double rate; // Instantaneous rate inside a window
};

/* Create another structure that contains coincidence data.
 * This has the number of photons on both PMTs plus total coinc length.*/
struct coinc_t { 
	// From Data
	unsigned long time;
	double realtime;
	// From generating UCN
	int pmt1;
	int pmt2;
	int prompt; // Number in the first window
	double length;
	double avg1;
	double avg2;
	double rate; // Instantaneous rate
	
	unsigned long index; // For loading coincidence PMT vectors (since they're saved)
};

struct counting_t {
	
	double dt;
	double N;

};

/* Create a Measurement Struct, which contains the values and errors 
 * of our varous run inputs. */
struct measurement {
	double val;
	double err;
	measurement operator/(const measurement& rhs) {
		return measurement {val/rhs.val, val/rhs.val*sqrt(pow(err/val, 2.0)+pow(rhs.err/rhs.val, 2.0))};
	}
	measurement operator*(const measurement& rhs) {
		return measurement {val*rhs.val, val*rhs.val*sqrt(pow(err/val, 2.0)+pow(rhs.err/rhs.val, 2.0))};
	}
	measurement operator+(const measurement& rhs) {
		return measurement {val+rhs.val, sqrt(err*err+rhs.err*rhs.err)};
	}
	measurement operator-(const measurement& rhs) {
		return measurement {val-rhs.val, sqrt(err*err+rhs.err*rhs.err)};
	}
};

/* Creating the RUN class. This contains the data from the run */
class Run {
  
  private:
	//Some of the files seem to have times assigned to the IO register that are twice as slow as they should be. Needs correction with ioMult (?)
	
	/* in Run.cpp */
	// Arguments for a given RUN (input when RUN loads)
	int coincMode = 1; //1 for fixed-window; 2 for moving-window
	int coincWindow; // Proximity of two gammas required for a coincidence finder to start
	int peSumWindow; // Proximity of gammas to continue the coincidence finder
	int promptWindow; // Proximity of gammas for a "prompt" window
	int peSum; // Number of PEs required for a coincidence
	int runNo; // Run number 
	
	int coinChan1 = 1; // 2-part coincidence channel
	int coinChan2 = 2;
	int mcsOffset = 0; // initial offset of MCS channel (MCS2 changes to +10 at 9020)
	
	unsigned long fastDT; // Hardcoded in -- coincidence "Fast Deadtime" (1.5* real deadtime)
	
	// Data loaded from ROOT by RUN
	char* fileName;
	TFile* dataFile;
	
	// Vectors generated by Run-makeHists, but initialized here
	std::vector<input_t> data;    //A vector holding all events in the file for RUN
	//std::vector<input_t> coinc;
	std::vector<coinc_t> coinc;
	//std::vector<std::vector<input_t> > pmtACoincHits;
	//std::vector<std::vector<input_t> > pmtBCoincHits;		
	std::vector<std::vector<input_t> > acCoincHits;
	
	// ROOT histograms of waveforms and photons spectra
	//TH1D phsA;
	//TH1D phsB;
	//TH1D pmt1SummedWaveform;
	//TH1D pmt2SummedWaveform;
			
	/* in Run-findcoincidence.cpp */
	void findCoincidenceFixed();
	void findCoincidenceMoving();
	void findCoincidenceFixedTele(); // Realistically should be able to set "PESum" to 0
	//void findCoincidenceFixedHT();
	//void findCoincidenceMovingHT();
	void findCoincidenceAC();
	void findAntiCoincidenceAC();
	void findCoincidenceNoPMT();
	//void findFixedHistos();
	//void findMovingHistos();
	
	/* in Run-readDataRoot.cpp */
	int numBits(uint32_t i);
	//void readDataRoot();
	void readDataRoot(const char* namecycle);
		
  public:
	
	/* in Run.cpp */
	//Run(int coincWindow, int peSumWindow, int peSum, const char* fName, int coincMode);
	//Run(int coincWindow, int peSumWindow, int peSum, int runNo, int coincMode, std::string runBody);
	Run();
	Run(int coincWindow, int peSumWindow, int peSum, int runNo, int coincMode, std::string runBody, const char* namecycle);
	//Run(int coincWindow, int peSumWindow, int peSum, std::vector<input_t> cts, int coincMode);
	~Run();
	
	bool exists();
	
	int getCoincMode();
	int getCoincWindow();
	int getPeSumWindow();
	int getPeSum();
	int getRunNo();
	
	int getCoinC1();
	int getCoinC2();
	int getMCSOff();
	
	
	std::vector<std::vector<input_t> > pmtACoincHits;
	std::vector<std::vector<input_t> > pmtBCoincHits;
	std::vector<std::vector<input_t> > pmtCCoincHits; // Sorted coinc hits
	//int coinChan1 = 1; // 2-part coincidence channel
	//int coinChan2 = 2;
	//int mcsOffset = 0; // initial offset of MCS channel (MCS2 changes to +10 at 9020)
	
	//TH1D getpmt1Waveform();
	//TH1D getpmt2Waveform();
	//TH1D getphsA();
	//TH1D getphsB();
	
	void setCoincMode(int mode);
	void setCoincWindow(int window);
	void setPeSum(int sum);
	void setPeSumWindow(int window);
	
	void setCoinC1(int chan);
	void setCoinC2(int chan);
	void setMCSOff(int off);
	
	/* in Run-getTagBitEvt.cpp */
	double getTagBitEvt(int mask, double offset, bool edge);
	
	/* in Run-makeHists.cpp */
	std::vector<input_t> getCounts(const std::function <input_t (input_t)>& expr, const std::function <bool (input_t)>& selection);
	//std::vector<input_t> getCoincCounts(const std::function <input_t (input_t)>& expr, const std::function <bool (input_t)>& selection);
	std::vector<coinc_t> getCoincCounts(const std::function <coinc_t (coinc_t)>& expr, const std::function <bool (coinc_t)>& selection);
	std::vector<size_t> getCoincIndices(const std::function<size_t (size_t)>& expr, const std::function <bool (coinc_t)>& selection);
	//std::vector<input_t> getCoincCountsAC(const std::function <input_t (input_t)>& expr, const std::function <bool (input_t)>& selection);
	std::vector<coinc_t> getCoincCountsAC(const std::function <coinc_t (coinc_t)>& expr, const std::function <bool (coinc_t)>& selection);
	//std::vector<double> getPhotonTracesVect(int pmt, const std::function <double (std::vector<input_t>)>& expr,const std::function <bool (std::vector<input_t>)>& selection);
	std::vector<counting_t> getCountingAvgs(int nvar, const std::function <coinc_t (coinc_t)>& expr,
						const std::function <bool (coinc_t)>& selection);
	//TH1D getCoincHist(const std::function <double (input_t)>& expr, const std::function <bool (input_t)>& selection);
	//TH1D getDeadtimeHist(double start, double end);
	//TH1D getHist(const std::function <double (input_t)>& expr, const std::function <bool (input_t)>& selection);
	//TH1D getHistIterator(
	//	const std::function <double (std::vector<input_t>::iterator, std::vector<input_t>::iterator, std::vector<input_t>::iterator)>& expr, 
	//	const std::function <bool (std::vector<input_t>::iterator, std::vector<input_t>::iterator, std::vector<input_t>::iterator)>& selection);
	//std::vector<double> allpmtAHits;
	//std::vector<double> allpmtBHits;
		
};

/* Pulling out the old commented code and putting it here for later:
 * //	double vUp;  //Holds the realtime of the vanadium going up for counting in the long run
//	double tdUp;
* //double clUp;
 * //	static bool before(input_t a, input_t b);
//	void findTDup();
//	void findVup();
//	void findCLup();
 * //	measurement getGVNorm();
 * /*void setTimeWindow(int low, int high);
	void setEnergyWindow(int low, int high);
	void setioMult(int mult);
	void setFileName(char *fname);*/