#include <vector>
#include <algorithm>
#include <string>
#include <math.h>
#include "stdio.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMySQLServer.h"
#include "TMySQLResult.h"
#include "TSQLRow.h"
#include "Run.hpp"

/*	------------------------------------------------------------------------------------------------
	Author: Nathan B. Callahan
	Editor: Frank M. Gonzalez
	
	This class handles communication with a MySQL database which holds the run numbers and timings.
	It then iterates over the list of runs with its various methods and returns aggregate results.
	
	There is 1 constructor: It takes the parameters to be loaded for each Run object.
	
	It also accepts an argument that queries the MySQL database to get the list of runs.
	
	The method getRuns queries the database and prepares the list of runs (runs) and a list of X
	values to plot against (xs).
	
	The method getMeasurements accepts a function analyzer as its argument. It iterates through the
	list of runs and returns a vector which is populated by pushing the return value of applying the
	analyzer function to each run in the list.
	
	The method getXs returns the vector of X values from the MySQL query.
	
	The method sumHistograms accepts a function summer as well as histogram binning information.
	It creates a new histogram with the given size and sums up all the histograms given by applying
	summer to each of the runs in the list. It returns the summed histogram.
	------------------------------------------------------------------------------------------------	*/

#pragma once

/* class DBHandler holds the information in the database, and the information
 * on the MySQLServer is saved in this class too. */
 
class DBHandler
{
	
	private:
	TMySQLServer* runLogserv;
	TMySQLServer* emsServ;
	char* query;
	std::vector<int> runs;
	std::vector<std::string> runBodies;
	std::vector<double> xs;
	int coincWindow;
	int peSumWindow;
	int peSum;
	int coincMode;
	void getRuns();
	
	public:
	DBHandler(const char* sqlQuery, int coincWindow, int peSumWindow, int peSum, int coincMode);
	~DBHandler();
	std::vector<measurement> getMeasurements(const std::function <measurement (Run*)>& analyzer);
	std::vector<double> getXs();
	TH1D sumHistograms(const std::function <TH1D (Run*)>& summer, int nbins, double low, double high);
	void foreach(const std::function <void (Run*)>& func);
	
};

/* Extra code that's being pulled out 
 PRIVATE:
 * //TFile* file;
 * //std::vector<measurement> results;
 * //std::function <measurement (Run*)> analyzer;
 * //std::function <TH1D (Run*)> summer;
 PUBLIC:
 * /*DBHandler(const char* sqlQuery, int coincWindow, int peSumWindow, int peSum, const std::function <measurement (Run*)>& analyzer);
	DBHandler(const char* sqlQuery, int coincWindow, int peSumWindow, int peSum, const std::function <TH1D (Run*)>& summer);
 *  //void printResults();
	*/
