#include "../inc/DBHandler.hpp"
#include "../inc/Run.hpp"

/*	------------------------------------------------------------------------------------------------
	Author: Nathan B. Callahan
	Editor: Frank M. Gonzalez
	
	This file contains the iterating methods which analyze the data.
	------------------------------------------------------------------------------------------------	*/


/* Accept a function which will return a vector<measurement> and will be evaluated for each run. 
 * We will return the vector<measurement> results. */
std::vector<measurement> DBHandler::getMeasurements(const std::function <measurement (Run*)>& analyzer) {
	std::vector<measurement> results; 
	
	// If the run list is empty, populate it.
	if(runs.empty()) { 
		this->getRuns();
	}
	
	// Loop through the runs. Begin by initializing loop vars
	auto it = runs.begin();
	char runName[256];
	auto bodiesIt = runBodies.begin();
	
	// Create the filename and the runobject. Loop through the total number of runs.
	for(it = runs.begin(); it < runs.end(); it++) {	
		char* fileLoc = std::getenv("RUNDIR");
		snprintf(runName, 256, "%s/processed_output_%%05d.root", fileLoc);
		sprintf(runName, runName, (*it));
		//sprintf(runName, "/media/frank/FreeAgentDrive/UCNtau/2017/processed_output/processed_output_%05d.root", (*it));
		printf("Opening Run %05d\n", (*it));
		Run run(this->coincWindow, this->peSumWindow, this->peSum, *it, coincMode, *bodiesIt);
		printf("Set coinc mode %d\n", coincMode);
		bodiesIt++;
		
		// skip any nonexistent runs 
		if(!run.exists()) {
			printf("Skipping Run %05d\n", (*it));
			continue;
		}
		
		// call the analyzer on our run, and call back the results
		measurement mes = analyzer(&run); 
		results.push_back(mes);
	}
	
	return results;
	
}

/* Accept a function which will return a TH1D (ROOT histogram) that will be 
 * evaluated for each run and summed. Returns the TH1D object. */
TH1D DBHandler::sumHistograms(const std::function <TH1D (Run*)>& summer, int nbins, double low, double high) {
	
	// If the run list is empty, populate it.
	if(runs.empty()) {
		this->getRuns();
	}
	
	// Loop through the runs. Begin by initializing loop vars
	auto it = runs.begin();
	char runName[256];
	int i;
	
	// Initialize our histogram. The DBHandler::sumHistograms object has 
	// some argument inputs defining the bins.
	TH1D summedHist("summedHist", "summedHist", nbins, low, high); 
	
	for(it = runs.begin(); it < runs.end(); it++) {
		char* fileLoc = std::getenv("RUNDIR");
		snprintf(runName, 256, "%s/processed_output_%%05d.root", fileLoc);
		sprintf(runName, runName, (*it));
		//sprintf(runName, "/media/frank/FreeAgentDrive/UCNtau/2016-2017/processed_output_%05d.root", (*it));
		printf("Opening Run %05d\n", (*it));
		Run run(this->coincWindow, this->peSumWindow, this->peSum, runName, coincMode);
		
		// Apply the (histogram) summer function to our runs. Loop through 
		// and sum all histograms.
		TH1D hist = summer(&run); 
		for(i = 0; i < nbins; i++) {
			summedHist.Fill(i, hist.GetBinContent(i));
		}
	}
	return summedHist;
}

/* Accept a function that will also open the DBHandler, without as much 
 * information as previously. */
void DBHandler::foreach(const std::function <void (Run*)>& func) {
	
	// If the run list is empty, populate it.
	if(runs.empty()) { 
		this->getRuns();
	}
	
	// Loop through the runs. Begin by initializing loop vars
	auto it = runs.begin();
	char runName[256];
	auto bodiesIt = runBodies.begin();
	
	// Loop through the runs to act on each run with the requisite 
	// predefined function.
	for(it = runs.begin(); it < runs.end(); it++) {
		char* fileLoc = std::getenv("RUNDIR");
		snprintf(runName, 256, "%s/processed_output_%%05d.root", fileLoc);
		sprintf(runName, runName, (*it));
		//sprintf(runName, "/media/frank/FreeAgentDrive/UCNtau/2017/processed_output/processed_output_%05d.root", (*it));
		printf("Opening Run %05d\n", (*it));
		
		// Create run Object
		Run run(this->coincWindow, this->peSumWindow, this->peSum, *it, coincMode, *bodiesIt);
		bodiesIt++;
		func(&run);
	}
}

/* Extra lines of code that've been commented out
 * //sprintf(runName, "/Volumes/SanDisk/2015-2016/raw_data/Run%05d.root", (*it)); //Create the filename
 * /*auto resultIt = results.begin();
	auto xsIt = xs.begin();
	
	printf("----Data----\n");
	
	for(resultIt = results.begin(); resultIt < results.end(); resultIt++) {
		printf("%f %f %f %f\n", (*xsIt), (*resultIt).val, 0.0, (*resultIt).err);
		xsIt++;
	}
	
	printf("----Data----\n");*/
//sprintf(runName, "/Volumes/SanDisk/2015-2016/raw_data/Run%05d.root", (*it));
		//sprintf(runName, "/Volumes/Seagate/raw_data/Oct2015-Feb2016/Run%05d.root", (*it));
