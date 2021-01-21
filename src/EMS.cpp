#include "../inc/EMS.hpp"

/*------------------------------------------------------------------------
   Author: Frank M. Gonzalez
   
	This is the code to load information from the EMS.
	
	We're creating an "EMS" object that contains information from the EMS
	part of our run trees.
	
	EMS Data structure: [bullet points are channels]
	Device 1200 -- Multiscan (Tau pressure, both OoM 1e-7)
	* 50 (North_CC)
	* 52 (South_CC)
	Device 9001 -- EPICS (Beam monitors)
	* 1 (LX current -- Somewhere between 8 and 10 when beam's on)
	* 2 (TBCM1 voltage -- should be around -15 when the accelerator's working)
	Device 1004 -- Flexrax (source pressure)
	* 1 (Source pressure)
	* 2 (PPM volume pressure)
	Device 3311 -- PMT Temps 
	* 0 (PMT1)
	* 1 (PMT2)
	
	
------------------------------------------------------------------------*/

//----------------------------------------------------------------------
EMS::EMS() {
	this->runNo = 0;
		
	// allocate memory and initialize file/trees
	fileNameEMS = new char[256];
	//dataFileEMS = NULL;
	//dataTree = NULL;
	// clean the trees 
	dataFileEMS = new TFile("", "read");
	
}
/* Load a run to create the pmt waveforms. Requires windows, sums, names, modes. */
EMS::EMS(int runNo, std::string runBody, const char* namecycle) {
	
	// save the input data into a tree (this) 
	this->runNo = runNo;
	
	// allocate memory and initialize file/trees
	fileNameEMS = new char[256];
	sprintf(fileNameEMS, runBody.c_str(), runNo);
	
	// clean the trees 
	//dataFileEMS = NULL;
	//clUp = 0.0;
	
	std::ifstream infile(fileNameEMS);
	if (infile.good()) {
		// load our data from file 
		dataFileEMS = new TFile(fileNameEMS, "read");
		//dataTree = NULL;
		//coincTree = NULL;
		
		// output of our summed waveforms, both in the pmts and in general.
		//pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
		//pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
		//phsA = TH1D("phsA", "phsA", 100, 0, 100);
		//phsB = TH1D("phsB", "phsB", 100, 0, 100);
		
		// check to read our data into manipulatable root tree
		this->readDataRoot(namecycle);
		if(data.empty()) {
			input_e blank;
			data.push_back(blank);
		}
		
		dataFileEMS->Close();
		//dataFileEMS = NULL;
	} //else {
	//	printf("Warning! Failed to load EMS data for Run %05d! \n", runNo);
	//}
	//fileName = NULL;
	//dataFile = NULL;
}

/* Destructor to clear Run data (saves memory)*/ 
EMS::~EMS() {	
	if(fileNameEMS != NULL) {	
		delete[] fileNameEMS;
	}
	if(dataFileEMS != NULL) {
		delete dataFileEMS;
	}
}


//----------------------------------------------------------------------
/* Initialization instructions. Stolen from RUN object. */
int EMS::numBits(uint32_t i)
{
     // Didn't come up with this, I assume it's magic. --NBC
     i = i - ((i >> 1) & 0x55555555);
     i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
     return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

/* Find the run number for our file */
int EMS::getRunNo() {
	return runNo;
}

/* Check to make sure the run actually loads normally */
bool EMS::exists() {
	if(dataFileEMS!=NULL){
		return(!dataFileEMS->IsZombie());
	}else{
		return(false);
	}
}

//----------------------------------------------------------------------
/* EMS data reading functions. We need to convert from ROOT to C++ */

/* Convert Root data into a C++ structure that we can actually use */
void EMS::readDataRoot(const char* namecycle) {
	// initialize variables and ROOT tree 
	int numEntries;
	input_e event;
	input_e eve0; // initial event.
	int i;
	TTree* rawData = NULL;
	
	// make sure we've initialized EMS properly
	if(this->exists() == false) {
		return;
	}
	
	// call the TList from the data file, to check the elements of our root tree 
	if (dataFileEMS==NULL){ return; }
	TList* list = dataFileEMS->GetListOfKeys();
	
	if(list == NULL) {
		return;
	}	
	if(list->Contains(namecycle)) {	
		// Load our data from file
		dataFileEMS->GetObject(namecycle, rawData);
		// check if we got a tree from file
		if(rawData != NULL) {
			// load leaves from raw data tree
			TLeaf* timelf = rawData->FindLeaf("time");
			TLeaf* devlf = rawData->FindLeaf("device");
			TLeaf* chlf = rawData->FindLeaf("channel");
			TLeaf* vallf = rawData->FindLeaf("val");
			// transfer the leaves to the C++ event
			timelf->SetAddress(&event.time);
			devlf->SetAddress(&event.dev);
			chlf->SetAddress(&event.ch);
			vallf->SetAddress(&event.val);
		// if no tree, return with empty vectors
		} else { 
			return;
		}

		// load the number of entries from the raw data and loop through
		// them. We can then find the events associated with each entry
		numEntries = rawData->GetEntries();
		for(i = 0; i < numEntries; i++) {
			rawData->GetEntry(i);
			data.push_back(event);
			
			//printf("%d,%d,%d,%e\n",event.time,event.ch,event.dev, event.val);
		}
		// sort the data, assuming we have data
		if(!data.empty()) {
			std::sort(data.begin(), data.end(), [](input_e x, input_e y)->bool{return (x.time < y.time);});
		}
	}
	return;
}

//----------------------------------------------------------------------
/* get the average on a dev/channel in a time interval*/ 
double EMS::getEMSAvg(int device, int channel, double cntStart, double cntEnd) {
	
	// Counting averages
	int cts = 0;
	double sum = 0.0;
	double avg;
	
	// Make sure we have data -- else return -1
	// Word of warning -- TBCM is negative (it's the only EMS event that is)
	if(data.empty()) {
		return -1;
	}
	
	// Load our event, loop through events and pick out appropriate ones.
	input_e event;
		
	for (auto stepIt = data.begin(); stepIt <= data.end(); stepIt++) {
		event = *stepIt;
		
		// filter device/channel
		if (event.dev == device && event.ch == channel) {
			// only pick events in appropriate time.
			// Note that the time from EMS is an integer. I cast here because I'll forget in functions.
			if (((double)event.time > cntStart) && ((double)event.time < cntEnd)) {
				sum += event.val;
				cts += 1;
			}
		}
	}
	// If we got counts in our range, return average. Else return -1.
	if (cts > 0) {
		avg = sum/(double)cts;
		return avg;
	} else {
		//printf("\n\nWarning! Device %d, channel %d doesn't have EMS data in this interval!\n",device,channel);
		return -1;
	}
}

double EMS::getEMSMax(int device, int channel, double cntStart, double cntEnd) {
	
	// Counting averages
	int cts = 0;
	double max = 0.0;
	
	// Make sure we have data -- else return -1
	// Word of warning -- TBCM is negative (it's the only EMS event that is)
	if(data.empty()) {
		return -1;
	}
	
	// Load our event, loop through events and pick out appropriate ones.
	input_e event;
		
	for (auto stepIt = data.begin(); stepIt <= data.end(); stepIt++) {
		event = *stepIt;
		
		// filter device/channel
		if (event.dev == device && event.ch == channel) {
			// only pick events in appropriate time.
			// Note that the time from EMS is an integer. I cast here because I'll forget in functions.
			if (((double)event.time > cntStart) && ((double)event.time < cntEnd)) {
				if (max < event.val) {
					max = event.val;
				}
				cts += 1;
			}
		}
	}
	// If we got counts in our range, return average. Else return -1.
	if (cts > 0) {
		return max;
	} else {
		//printf("\n\nWarning! Device %d, channel %d doesn't have EMS data in this interval!\n",device,channel);
		return -1;
	}
}

double EMS::getEMSMin(int device, int channel, double cntStart, double cntEnd) {
	// Counting averages
	int cts = 0;
	double min = 1.0; // Assume peak is < 1
	
	// Make sure we have data -- else return -1
	// Word of warning -- TBCM is negative (it's the only EMS event that is)
	if(data.empty()) {
		return -1;
	}
	
	// Load our event, loop through events and pick out appropriate ones.
	input_e event;
			
	for (auto stepIt = data.begin(); stepIt <= data.end(); stepIt++) {
		event = *stepIt;
		
		// filter device/channel
		if (event.dev == device && event.ch == channel) {
			// only pick events in appropriate time.
			// Note that the time from EMS is an integer. I cast here because I'll forget in functions.
			if (((double)event.time > cntStart) && ((double)event.time < cntEnd)) {
				if(event.val < min) { 
					min = event.val;
				}
				cts += 1;
			}
		}
	}
	// If we got counts in our range, return average. Else return -1.
	if (cts > 0) {
		return min;
	} else {
		return -1;
	}
}


/* Find a vector of counts from our data set */
std::vector<input_e> EMS::getCounts(const std::function <input_e (input_e)>& expr,
									const std::function <bool (input_e)>& selection)
{
	// load our data vectors
	std::vector<input_e> filtered = {};
	std::vector<input_e> transformed = {};
	
	// If we haven't loaded our root tree, return.
	// Probably will give us a segfault since it returns a null vector.
	if(!data.empty()) {
		// copy and transform our data sets to find the total amount of counts
		std::copy_if(data.begin(), data.end(), std::back_inserter(filtered), selection);
		std::transform(filtered.begin(), filtered.end(), std::back_inserter(transformed), expr);
	} //else {
	//	printf("Ah fuck I can't believe you've done this\n");
	//}
	return transformed;
}
//runMCS2(coincWindow, peSumWindow, peSum, runNo, coincMode, fName, "tmcs_1");

/* Here we begin creating the output structures. These short functions 
 * are called in the main Run() function and each generate a pmt */
 
 /* Histograms that contain the waveform from each PMT. */
/*TH1D EMS::getpmt1Waveform() {
	if(coinc.empty()) {
		if(coincMode == 1) {
		//	this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
		//	this->findcoincidenceMoving();
		}
	}
	return pmt1SummedWaveform;
}
TH1D EMS::getpmt2Waveform() {
	if(coinc.empty()) {
		if(coincMode == 1) {
		//	this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
		//	this->findcoincidenceMoving();
		}
	}
	return pmt2SummedWaveform;
}

/* Histograms that contain the waveform for each position. */
/*TH1D EMS::getphsA() {
	if(coinc.empty()) {
		if(coincMode == 1) {
		//	this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
		//	this->findcoincidenceMoving();
		}
	}
	return phsA;
}
TH1D EMS::getphsB() {
	if(coinc.empty()) {
		if(coincMode == 1) {
		//	this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
		//	this->findcoincidenceMoving();
		}
	}
	return phsB;
}



//----------------------------------------------------------------------
/* Choose some parameters of the run cycle. We reset the counters each
 * time and then set the value we want. */
/*void EMS::setCoincWindow(int window) {
	coincWindow = window;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmt1SummedWaveform.Reset();
	pmt2SummedWaveform.Reset();
}
void EMS::setPeSumWindow(int window) {
	peSumWindow = window;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmt1SummedWaveform.Reset();
	pmt2SummedWaveform.Reset();
}
void EMS::setPeSum(int sum) {
	peSum = sum;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmt1SummedWaveform.Reset();
	pmt2SummedWaveform.Reset();
}
void EMS::setCoincMode(int mode) {
	coincMode = mode;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmt1SummedWaveform.Reset();
	pmt2SummedWaveform.Reset();
}

//----------------------------------------------------------------------
/* Readout for our different data parameters. */
/*int EMS::getCoincWindow() {
	return coincWindow;
}
int EMS::getPeSumWindow() {
	return peSumWindow;
}
int EMS::getPeSum() {
	return peSum;
}
int EMS::getCoincMode() {
	return coincMode;
}


/* Load a run to create the pmt waveforms. Here we're creating an arbitrary
 * file name, loading it from file*/
/*EMS::EMS(int coincWindow, int peSumWindow, int peSum, int runNo, int coincMode, std::string runBody) {
	
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
	clUp = 0.0;
	
	// load our data from file
	dataFile = new TFile(fileName, "read");
	dataTree = NULL;
	coincTree = NULL;
	
	// output of our summed waveforms, both in the pmts and in general.
	pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	phsA = TH1D("phsA", "phsA", 100, 0, 100);
	phsB = TH1D("phsB", "phsB", 100, 0, 100);
}
*/
/*EMS::EMS(int coincWindow, int peSumWindow, int peSum, const char* fName, int coincMode) {
	
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
	clUp = 0.0;
	
	// load our data from file 
	dataFile = new TFile(fName, "read");
	fileName = strdup(fName);
	dataTree = NULL;
	coincTree = NULL;
	
	// output of our summed waveforms
	pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
}*/

/* Load a run to create our pmt and our general waveforms. Now we also 
 * plan on loading this into ROOT at the end. */
/*EMS::EMS(int coincWindow, int peSumWindow, int peSum, int runNo, int coincMode, std::string runBody, const char* namecycle) {
	
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
	clUp = 0.0;
	
	std::ifstream infile(fileName);
	if (infile.good()) {
		// load our data from file 
		dataFile = new TFile(fileName, "read");
		dataTree = NULL;
		coincTree = NULL;
		
		// output of our summed waveforms, both in the pmts and in general.
		pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
		pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
		phsA = TH1D("phsA", "phsA", 100, 0, 100);
		phsB = TH1D("phsB", "phsB", 100, 0, 100);
		
		// check to read our data into manipulatable root tree
		//this->readDataRoot(namecycle);
		if(data.empty()) {
			input_t blank;
			data.push_back(blank);
		}
	} else {
		printf("Warning! Run %05d does not exist! \n", runNo);
	}
}

/* */
/*EMS::EMS(int coincWindow, int peSumWindow, int peSum, std::vector<input_t> cts, int coincMode) {
	
	// save the input data into a tree (this)
	this->coincWindow = coincWindow;
	this->peSumWindow = peSumWindow;
	this->peSum = peSum;
	this->runNo = 5;
	this->coincMode = coincMode;
	
	// clear the ROOT files 
	clUp = 0.0;
	fileName = NULL;
	dataFile = NULL;
	dataTree = NULL;
	coincTree = NULL;
	
	// load our input data
	data = cts;
	
	// output of our summed waveforms
	pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	phsA = TH1D("phsA", "phsA", 100, 0, 100);
	phsB = TH1D("phsB", "phsB", 100, 0, 100);
}

*/
