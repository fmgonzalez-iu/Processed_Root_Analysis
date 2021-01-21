#include "../inc/Run.hpp"

/*------------------------------------------------------------------------
   Author: Nathan B. Callahan
   Editor: Frank M. Gonzalez
   
   This is the code for all the "Run" functions. We use these to produce
   histograms representing the coincidence events of our analysis sims.
    
   There are 2 MCS boxes presently. Channel list:
   MCS1:
   * 1) Dagger PMT1
   * 2) Dagger PMT2
   * 3) Old GV
   * 4) Bare GV (RH top)
   * 5) SP
   MCS2: [Note that runs before 9020 are offset by 10]
   * 11) AC
   * 12) Downstream
   * 13) Foil (RH bot)
   * 14) Dagger PMT1 (high threshold)
   * 15) Dagger PMT2 (high threshold) 
   
------------------------------------------------------------------------*/
Run::Run() { // Blank Initiator
	// save the input data into this run object
	this->coincWindow = 0;
	this->peSumWindow = 0;
	this->promptWindow = 0;
	this->peSum = 0;
	this->runNo = 0;
	this->coincMode = 0;
	this->coinChan1 = 0;
	this->coinChan2 = 0;
	this->mcsOffset = 0;
		
	// allocate memory and initialize file/trees
	fileName = new char[256];
	
	// clean the trees 
	//dataFile = new TFile("", "read");
	
}


/* Load a run to create our pmt and our general waveforms. Now we also 
 * plan on loading this into ROOT at the end. */
Run::Run(int coincWindow, int peSumWindow, int peSum, int runNo, int coincMode, std::string runBody, const char* namecycle) {
	
	// save the input data into this run object
	this->coincWindow = coincWindow;
	this->peSumWindow = peSumWindow;
	this->peSum = peSum;
	this->runNo = runNo;
	this->coincMode = coincMode;
	// TODO: Update promptWindow to read input from file
	this->promptWindow = peSumWindow; // For now....
	
	// set initial coincidence channels based on mode
	if ((coincMode % 2) == 1) { // Odd for low threshold
		this->coinChan1 = 1;
		this->coinChan2 = 2;
	} else if ((coincMode % 2) == 0) { // Even for high threshold
		this->coinChan1 = 4;
		this->coinChan2 = 5;
	} else { // technically only a negative coincMode can do this
		printf("Warning! coincMode %d is invalid! Please modify the Run.cpp file!",coincMode);
		this->coinChan1 = 0;
		this->coinChan2 = 0;
	}
	
	if (coincMode == 9 || coincMode == 10) { // Removing fast coincidences
		this->fastDT = 30;
	} else { this->fastDT = 0; } // Set to 0 if we're doing "remove fast"
	
	// offset is set in the MAIN function since it's a year
	this->mcsOffset = 0;
		
	// allocate memory and initialize file/trees
	fileName = new char[256];
	sprintf(fileName, runBody.c_str(), runNo);
	
	// clean the trees 
	std::ifstream infile(fileName);
	if (infile.good()) {
		// load our data from file 
		dataFile = new TFile(fileName, "read");
				
		// check to read our data into manipulatable root tree
		this->readDataRoot(namecycle);
		if(data.empty()) {
			input_t blank;
			data.push_back(blank);
			printf("Warning! Run %05d, tree %s is empty! \n", runNo, namecycle);
		}
				
		// Number of photons in each PMT
		//phsA = TH1D("phsA", "phsA", 100, 0, 100);
		//phsB = TH1D("phsB", "phsB", 100, 0, 100);
		// Arrival time of photons
		//pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
		//pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
		dataFile->Close();
		
	} else {
		printf("Warning! Run %05d does not exist! \n", runNo);
	}
}

/* Destructor to clear Run data (saves memory)*/ 
Run::~Run() {
	if(fileName != NULL) {
		delete[] fileName;
	}
	if(dataFile != NULL) {
		delete dataFile;
	}
}

/* Find the run number for our file */
int Run::getRunNo() {
	return runNo;
}

//----------------------------------------------------------------------
/* Choose some parameters of the run cycle. We reset the counters each
 * time and then set the value we want. */
void Run::setCoincWindow(int window) {
	coincWindow = window;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmtCCoincHits.clear();
	acCoincHits.clear();
	//pmt1SummedWaveform.Reset();
	//pmt2SummedWaveform.Reset();
}
void Run::setPeSumWindow(int window) {
	peSumWindow = window;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmtCCoincHits.clear();
	acCoincHits.clear();
	//pmt1SummedWaveform.Reset();
	//pmt2SummedWaveform.Reset();
}
void Run::setPeSum(int sum) {
	peSum = sum;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmtCCoincHits.clear();
	acCoincHits.clear();
	//pmt1SummedWaveform.Reset();
	//pmt2SummedWaveform.Reset();
}
void Run::setCoincMode(int mode) {
	coincMode = mode;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmtCCoincHits.clear();
	acCoincHits.clear();
	//pmt1SummedWaveform.Reset();
	//pmt2SummedWaveform.Reset();
}

void Run::setCoinC1(int chan) {
	coinChan1 = chan;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmtCCoincHits.clear();
	acCoincHits.clear();
	//pmt1SummedWaveform.Reset();
	//pmt2SummedWaveform.Reset();
}

void Run::setCoinC2(int chan) {
	coinChan2 = chan;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmtCCoincHits.clear();
	acCoincHits.clear();
	//pmt1SummedWaveform.Reset();
	//pmt2SummedWaveform.Reset();
}

void Run::setMCSOff(int off) {
	mcsOffset = off;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmtCCoincHits.clear();
	acCoincHits.clear();
	//pmt1SummedWaveform.Reset();
	//pmt2SummedWaveform.Reset();
}
//----------------------------------------------------------------------
/* Readout for our different data parameters. */
int Run::getCoincWindow() {
	return coincWindow;
}
int Run::getPeSumWindow() {
	return peSumWindow;
}
int Run::getPeSum() {
	return peSum;
}
int Run::getCoincMode() {
	return coincMode;
}

int Run::getCoinC1() {
	return coinChan1;
}
int Run::getCoinC2() {
	return coinChan2;
}
int Run::getMCSOff() { 
	return mcsOffset;
}


/* Check to make sure the run actually loads normally */
bool Run::exists() {
	if(dataFile!=NULL){
		return(!dataFile->IsZombie());
	}else{
		return(false);
	}
}

//----------------------------------------------------------------------
/* Extra code for cleanliness.
 * 
 * 
 * /* Load a run to create the pmt waveforms. Requires windows, sums, names, modes. */
 //----------------------------------------------------------------------
/* Here we begin creating the output structures. These short functions 
 * are called in the main Run() function and each generate a histo*/
 
 /* Histograms that contain the number of photons for each PMT. */
/*TH1D Run::getphsA() {
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findCoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findCoincidenceMoving();
		}
		else if(coincMode == 3) {
			this->findCoincidenceFixedHT();
		}
		else if(coincMode == 4) {
			this->findCoincidenceMovingHT();
		}
	}
	return phsA;
}
TH1D Run::getphsB() {
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findCoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findCoincidenceMoving();
		}
		else if(coincMode == 3) {
			this->findCoincidenceFixedHT();
		}
		else if(coincMode == 4) {
			this->findCoincidenceMovingHT();
		}
	}
	return phsB;
}*/
 
 /* Histograms that contain the waveform from each PMT. */
/*TH1D Run::getpmt1Waveform() {
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findCoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findCoincidenceMoving();
		}
		else if(coincMode == 3) {
			this->findCoincidenceFixedHT();
		}
		else if(coincMode == 4) {
			this->findCoincidenceMovingHT();
		}
	}
	return pmt1SummedWaveform;
}
TH1D Run::getpmt2Waveform() {
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findCoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findCoincidenceMoving();
		}
		else if(coincMode == 3) {
			this->findCoincidenceFixedHT();
		}
		else if(coincMode == 4) {
			this->findCoincidenceMovingHT();
		}
	}
	return pmt2SummedWaveform;
}*/
/*Run::Run(int coincWindow, int peSumWindow, int peSum, const char* fName, int coincMode) {
	
	// save the input data into a tree (this)
	this->coincWindow = coincWindow;
	this->peSumWindow = peSumWindow;
	this->peSum = peSum;
	this->coincMode = coincMode;
	
	// reset the run number 
	runNo = -1;
	
	// load/clear the new file and required trees.
	strcpy(fileName, fName);
	dataFile = NULL;
	
	// load our data from file 
	dataFile = new TFile(fName, "read");
	fileName = strdup(fName);
	
	// Number of photons in each PMT
	phsA = TH1D("phsA", "phsA", 100, 0, 100);
	phsB = TH1D("phsB", "phsB", 100, 0, 100);
	// Arrival time of photons
	pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
}*/

/* Load a run to create the pmt waveforms. Here we're creating an arbitrary
 * file name, loading it from file*/
/*Run::Run(int coincWindow, int peSumWindow, int peSum, int runNo, int coincMode, std::string runBody) {
	
	// save the input data into a tree (this)
	this->coincWindow = coincWindow;
	this->peSumWindow = peSumWindow;
	this->peSum = peSum;
	this->runNo = runNo;
	this->coincMode = coincMode;
	
	// allocate memory and initialize file/trees
	fileName = new char[256];
	sprintf(fileName, runBody.c_str(), runNo);
	
	// clean the trees 
	dataFile = NULL;
		
	// load our data from file
	dataFile = new TFile(fileName, "read");
	
	// Number of photons in each PMT
	phsA = TH1D("phsA", "phsA", 100, 0, 100);
	phsB = TH1D("phsB", "phsB", 100, 0, 100);
	// Arrival time of photons
	pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	
}*/

//fileName = fName;
/* //gvNorm = (measurement){0.0, 0.0};
	//vUp = 0.0;
	//tdUp = 0.0;
	
	//sprintf(fileName, "/Volumes/Seagate/raw_data/Oct2015-Feb2016/Run%05d.root", runNo);
	//fileName = fName;
	//gvNorm = (measurement){0.0, 0.0};
	//vUp = 0.0;
	//tdUp = 0.0;
	
	//sprintf(fileName, "/Volumes/Seagate/raw_data/Oct2015-Feb2016/Run%05d.root", runNo);
	//fileName = fName;
	//gvNorm = (measurement){0.0, 0.0};
	//vUp = 0.0;
	//tdUp = 0.0;
	*/



/* */
/*Run::Run(int coincWindow, int peSumWindow, int peSum, std::vector<input_t> cts, int coincMode) {
	
	// save the input data into a tree (this)
	this->coincWindow = coincWindow;
	this->peSumWindow = peSumWindow;
	this->peSum = peSum;
	this->runNo = 5;
	this->coincMode = coincMode;
	
	// clear the ROOT files 
	fileName = NULL;
	dataFile = NULL;

	// load our input data
	data = cts;
	
	// Number of photons in each PMT
	phsA = TH1D("phsA", "phsA", 100, 0, 100);
	phsB = TH1D("phsB", "phsB", 100, 0, 100);
	// Arrival time of photons
	pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
}*/
