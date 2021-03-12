//#include "inc/DBHandler.hpp"
#include "inc/Run.hpp"
#include "inc/EMS.hpp"
#include "inc/Pileup.hpp"
#include "inc/Functions.hpp"
#include "inc/ExpFill.hpp"
#include "inc/GlobalFiles.hpp"

//#include "TMySQLServer.h"
//#include "TMySQLResult.h"
//#include "TMySQLRow.h"
#include "TPad.h"
#include "TF1.h"
#include <iostream>
#include <string>
#include <sstream>

#define NANOSECOND .000000001

/* ---------------------------------------------------------------------
 * Author: Frank M. Gonzalez
 * 
 * This is the main function for our UCNtau analyzer. We initialize some
 * variables called from our run files and then read through our analyses.
 * 
 * This analyzer has a lot of functionality, and it's relatively modular.
 * For questions about this ask frangonz@indiana.edu
 * 
 * -------------------------------------------------------------------*/

// Loaded coincidence parameters
struct params_in {
	 
	int coincWindow;
	int peSumWindow;
	int peSum;
	int coincMode;
	int jobid; 
};

std::vector<int> load_run_list(std::string query) {
	std::istringstream iss(query); // convert argv[1] into stringstream
	std::string token;
	
	//int nRuns = 0; // Run counter
	std::vector<int> rList;
	while(std::getline(iss, token, ',')) { // This should load a file
		int runNo = atoi(token.c_str());
		if (runNo > 0) { // If we just entered a run or list of runs
			rList.push_back(runNo);
		} else {
			std::ifstream run_file(query.c_str()); // Treats the query as a file
			if (run_file.is_open()) { // Check if the run_file works
				printf("Able to load file:\n");
				std::string line;
				while(std::getline(run_file,line,',')) {
					rList.push_back(std::stoi(line));
					printf("%05d,",std::stoi(line));
				}
			}
		}
	}
	
	return rList;
};

char* getFileName(int run_num) {
	
	// Trying to get 2017 and 2018 to both work.
	const char* fileLoc17 = std::getenv("RUNDIR_2017");
	const char* fileLoc18 = std::getenv("RUNDIR_2018");
	const char* fileLoc19 = std::getenv("RUNDIR_2019");
	
	char* fName = new char[256];
	
	// Also hard coding in run 3412 as the beginning of 2017.
	// 9519-9538 are in both directories BUT only in 2018 runlogdump.
	//if ((3412 <= run_num) && (run_num < 9519)) {
	if ((3412 <= run_num) && (run_num < 9927)) {
		snprintf(fName, 256, "%s/processed_output_%%05d.root", fileLoc17);
	//} else if ((9519 <= run_num) && (run_num < 14732)) {
	} else if ((9927 <= run_num) && (run_num < 14732)) {
		snprintf(fName, 256, "%s/processed_output_%%05d.root", fileLoc18);
	} else if((14732 <= run_num) && (run_num < 17999)) {
		snprintf(fName, 256, "%s/processed_output_%%05d.root", fileLoc19);
	} else {
		printf("Warning! Run %05d is not from calendar year 2017 or 2018 !\n", run_num);
		return fName;
	}
	
	std::string runBody = fName;
	sprintf(fName, runBody.c_str(), run_num);
	
	return fName;
}
	
void functionality_extract(std::vector<int> runList, params_in params, int setting) {
	// Function Extraction
	const char* fillLoc  = std::getenv("FILL_LOC");  // MAD fill fit time const.
	const char* traceLoc = std::getenv("TRACE_LOC"); // Exp fill time const.
	const char* saveLoc  = std::getenv("SAVE_LOC");  // Save location for ROOT drawings
	
	char* fName;
	for (auto rIt = runList.begin(); rIt < runList.end(); rIt++) {
	
		int runNo = *rIt;
		fName = getFileName(runNo);
		std::ifstream infile(fName);
		if (infile.good()) {
			Run runMCS1(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_0");
			Run runMCS2(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_1");
			EMS runEMS(runNo, fName, "tems");
			
			if (runNo >= 9020) { runMCS2.setMCSOff(10);	}
			switch(setting) {
				case 1: {normNByDip(&runMCS1, &runMCS2, &runEMS, traceLoc);     } break;
				case 2: {normNByDipSing(&runMCS1, &runMCS2, &runEMS, traceLoc); } break;
				case 3: {normNByDipCoinc(&runMCS1, &runMCS2, &runEMS, traceLoc);} break;
				case 4: {extractObservables(&runMCS1, &runMCS2, traceLoc);} break;
				case 5: {activeCleaner(&runMCS1,&runMCS2,&runEMS,traceLoc);} break;
				case 6: {fitBkgUnload(&runMCS1,&runMCS2,&runEMS);} break;
				case 7: {normNByDipDet(&runMCS1,&runMCS2,&runEMS,traceLoc);} break;
				case 8: {holdingPressure(&runMCS1,&runMCS2,&runEMS);} break;
				default:{
					printf("Invalid setting %d for function 1!\n",setting);
				} 
			}
		}
		fflush(stdout);
	}
}

void save1DRootCSV(char* saveNameTmp, TH1D hist, int nBins) {
	// Convert a 1D Root Histogram to a csv
	 
	char saveName[256]; // First, copy the 256 character pointer to the saveName
	strncpy(saveName,saveNameTmp,256);
	FILE* outfile = fopen(saveName,"w");
	for (int ii = 0; ii < (nBins - 1); ii++) { // And save
		if ((hist.GetBinContent(ii) >= 0) && (hist.GetBinContent(ii+1) >= 0)) { // Make sure it's positive
			fprintf(outfile,"%d,%d\n",ii,(int)(hist.GetBinContent(ii)));
		} else {
			fprintf(outfile,"%d,%d\n",ii,0);
		}
	}
	fclose(outfile);
}

void save1DRootCSV_float(char* saveNameTmp, TH1D hist, int nBins) {
	// Convert a 1D Root Histogram to a csv
	 
	char saveName[256]; // First, copy the 256 character pointer to the saveName
	strncpy(saveName,saveNameTmp,256);
	FILE* outfile = fopen(saveName,"w");
	for (int ii = 0; ii < (nBins - 1); ii++) { // And save
		fprintf(outfile,"%d,%f\n",ii,hist.GetBinContent(ii));
	}
	fclose(outfile);
}

void save3DRootCSV(char* saveNameTmp, TH1D pmt1, TH1D pmt2, TH1D coinc, int nBins) {
	// Convert 3 Root Histogram to a csv (pmt 1/2/C)
	
	char saveName[256];  // First, copy the 256 character pointer to the saveName
	strncpy(saveName, saveNameTmp, 256);
	FILE* outfile = fopen(saveName, "w");
	for (int ii = 0; ii < nBins; ii++) {
		fprintf(outfile,"\n%d",ii);
		fprintf(outfile,",%f",pmt1.GetBinContent(ii));
		fprintf(outfile,",%f",pmt2.GetBinContent(ii));
		fprintf(outfile,",%f",coinc.GetBinContent(ii));
	}
	fclose(outfile);
}

void functionality_summing(std::vector<int> runList, params_in params, int setting) {
// For summing up a bunch of individual PMT hits.
// We're using ROOT histograms for these, so that's what we'll want to generate first.

	const char* fillLoc  = std::getenv("FILL_LOC");  // MAD fill fit time const.
	const char* traceLoc = std::getenv("TRACE_LOC"); // Exp fill time const.
	const char* saveLoc  = std::getenv("SAVE_LOC");  // Save location for ROOT drawings
	
	if (saveLoc == NULL) {
		printf("SAVE_LOC not set! Defaulting to './' !\n");
		saveLoc = "./";
	}
	
	int nBins;
	int start;
	int end;
	int nBins2;
	int start2;
	int end2;
	// Dip summing -- setting (a lot...)
	//TH1D dipSumS("summedDipS", "summedDipS", 3100, 0, 310);
	
	// Sums of long holds -- Setting 7
	switch (setting) {
		case 1: { nBins=3100; start=0; end=310;  } break; // Unload, short
		case 2: { nBins=3100; start=0; end=310;  } break; // Unload, long
		case 3: { nBins=1;    start=0; end=1;    } break; // FitFillFree, doesn't need hists.
		case 4: { nBins=3100; start=0; end=310;  } break; // Fast Counter
		case 5: { nBins=1;    start=0; end=1;    } break; // FitFill, doesn't need hists
		case 6: { nBins=2200; start=0; end=2200; } break; // Active Cleaner Summing
		case 7: { nBins=1550; start=0; end=1550; } break; // Long Holding Times
		case 8: { nBins=1000; start=0; end=1000; } break; // Just putting this in to not brak thing
		case 9: { nBins=3100; start=0; end=310;  } break; // Summing Peak 1
		case 10:{ // Monitor Summer
			nBins=2200; start=0; end=2200; } break; // Put this in for full run monitors
			if (params.jobid >= 9000) { // This isn't the actual break, but for the memes.
				nBins=580;  start=0; end=580; // For 20s holds
			} else {  
				nBins=730;  start=0; end=730; } } break; 
		case 11:{ nBins=2200; start=0; end=2200; } break; // Full Run Summer (dagger)
		case 12:{ nBins=2600; start=0; end=260; } break; // Technically coincidence structure.
	}
	
	// There's 5 unique 1D histograms at any given time.
	TH1D hist1("summed1","summed1",nBins,start,end);
	TH1D hist2("summed2","summed2",nBins,start,end);
	TH1D hist3("summed3","summed3",nBins,start,end);
	TH1D hist4("summed4","summed4",nBins,start,end);
	TH1D hist5("summed5","summed5",nBins,start,end);
	TH1D hist6("summed6","summed6",nBins,start,end);
	char* fName;
	
	for (auto rIt = runList.begin(); rIt < runList.end(); rIt++) {
	
		int runNo = *rIt;
		fName = getFileName(runNo);
		std::ifstream infile(fName);
		if (infile.good()) {
			Run runMCS1(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_0");
			Run runMCS2(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_1");
			EMS runEMS(runNo, fName, "tems");
			
			if (runNo >= 9020) { runMCS2.setMCSOff(10);	}
			
			switch(setting) { 
				case 1: {dipSummerS(&runMCS1, &runMCS2, &hist1, 20.0); } break;
				case 2: {dipSummerL(&runMCS1, &runMCS2, &hist1); } break;
				case 3: {
					for (int ii = 1; ii < 11; ii++) {
						fitFillFree(&runMCS1,&runMCS2,ii);	
					}} break;
				//case 4: {fitLastDip(&runMCS1, &runMCS2, &dipSumS);} break;
				case 4: {writeFastEvts(&runMCS1,&runMCS2,&hist1);
						 writeFastBkgs(&runMCS1,&runMCS2,&hist2);} break;
				case 5: {
					for (int ii = 1; ii < 11; ii++) {
						fitFill(&runMCS1,&runMCS2,fillLoc,ii);	
					}} break;
				case 6: {acSummerCoinc(&runMCS1, &runMCS2, &hist1, 20.0);} break;
				case 7: {longHoldSummer(&runMCS1, &runMCS2, &hist1, &hist2, &hist3);} break;
				case 8: {bkgSummer(&runMCS1, &runMCS2, &hist1,&hist2,&hist3);} break;
				case 9: {peak1SummerS(&runMCS1, &runMCS2,&hist1,&hist2,true);
						 peak1SummerS(&runMCS1, &runMCS2,&hist3,&hist4,false);} break;
						 //peak1SummerL(&runMCS1, &runMCS2,&hist3,&hist4);} break;
				case 10: { bool multiHold = true;
						  if (params.jobid > 1000) { 	multiHold = false;}
						  monSummer(&runMCS1,&runMCS2,nBins,3,&hist1,multiHold); // GV
						  monSummer(&runMCS1,&runMCS2,nBins,4,&hist2,multiHold); // Ba/RHAC
						  monSummer(&runMCS1,&runMCS2,nBins,5,&hist3,multiHold); // SP
						  monSummer(&runMCS1,&runMCS2,nBins,7,&hist4,multiHold); // DS
						  monSummer(&runMCS1,&runMCS2,nBins,8,&hist5,multiHold); } break; // Foil/RH
				case 11: {totalSummer(&runMCS1, &runMCS2, &hist1, &hist2, &hist3, 20.0);}break;
				case 12: {numPhotonsByTime(&runMCS1, &runMCS2, &hist1, &hist2, &hist3, 20.);
						  numPhotonsByTime(&runMCS1, &runMCS2, &hist4, &hist5, &hist6, 1550.); } break;
				default: {
					printf("Invalid setting %d for function 2!\n",setting);
				} 
			}
		}
		fflush(stdout);
	}
	
	// Save histograms
	// Initialize saving names
	char dirName[256];
	snprintf(dirName, 256, "%s/%%s", saveLoc);
	printf("Saving Histograms and .CSVs in folder %s",dirName);
	char rootName[256];
	char tmpName[256];
	char csvName[256];
	switch(setting) {
		case 1: {
			snprintf(tmpName,256,dirName,"summedDipShort.root%d"); // Root
			snprintf(rootName,256,tmpName,params.jobid);
			hist1.SaveAs(rootName);
			
			snprintf(tmpName,256,dirName,"summedDip.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			save1DRootCSV(csvName,hist1,nBins);
			} break;
		case 2: {
			snprintf(tmpName,256,dirName,"summedDipLong.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist1.SaveAs(rootName);
			} break;
		case 3: { } break; // No saved for case 3
		case 4: {
			// Create character names
			if (hist1.GetBinContent(0) > 0) { // Saving Unload
				//printf("Saving Unload\n");
				snprintf(tmpName,256,dirName,"summedRemove.root%d");
				snprintf(rootName,256,tmpName,params.jobid);
				hist1.SaveAs(rootName);
				
				snprintf(tmpName,256,dirName,"oneRemove.csv%d");
				snprintf(csvName,256,tmpName,params.jobid);
				save1DRootCSV(csvName,hist1,nBins);
			}
			if (hist2.GetBinContent(0) > 0) { // Saving Background
				//printf("Saving Background\n");
				snprintf(tmpName,256,dirName,"summedRemoveBkg.root%d");
				snprintf(rootName,256,tmpName,params.jobid);
				hist2.SaveAs(rootName);
				
				snprintf(tmpName,256,dirName,"oneRemoveBkg.csv%d");
				snprintf(csvName,256,tmpName,params.jobid);
				save1DRootCSV(csvName, hist2, nBins);
			}
			} break;
		case 5: {	} break;
		case 6: {
			snprintf(tmpName,256,dirName,"acSumAll");
			snprintf(rootName,256,tmpName,params.jobid);
			hist1.SaveAs(rootName);
			} break;
		case 7: {
			// Save the root files
			snprintf(tmpName,256,dirName,"longHold1.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist1.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"longHold2.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist2.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"longHoldC.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist3.SaveAs(rootName);
			
			// Save an output .csv file
			snprintf(tmpName,256,dirName,"holdingBkg.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			save3DRootCSV(csvName, hist1, hist2, hist3, nBins);
			
			} break;
		case 8: {
			snprintf(tmpName,256,dirName,"bkgSum1.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist1.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"bkgSum2.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist2.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"bkgSumC.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist3.SaveAs(rootName);
			// Save an output .csv file
			snprintf(tmpName,256,dirName,"bkgByPMT.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			save3DRootCSV(csvName,hist1,hist2,hist3,nBins);
			} break;
		case 9: {
			snprintf(tmpName,256,dirName,"dipSumACS.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist2.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"dipSumACS.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			save1DRootCSV_float(csvName, hist2, nBins);
			
			snprintf(tmpName,256,dirName,"dipSumS.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist1.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"dipSumS.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			save1DRootCSV_float(csvName, hist1, nBins);
			
			snprintf(tmpName,256,dirName,"dipSumACL.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist4.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"dipSumACL.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			save1DRootCSV_float(csvName, hist4, nBins);
			
			snprintf(tmpName,256,dirName,"dipSumL.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist3.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"dipSumL.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			save1DRootCSV_float(csvName, hist3, nBins);
			
			}break;
		case 10: {
			// Hard code this stuff	
			snprintf(tmpName,256,dirName,"monSumGV_%05d.root");
			snprintf(rootName,256,tmpName,params.jobid);
			hist1.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"monSumGV_%05d.csv");
			snprintf(csvName,256,tmpName,params.jobid);
			save1DRootCSV_float(csvName, hist1, nBins);
			
			snprintf(tmpName,256,dirName,"monSumRHAC_%05d.root");
			snprintf(rootName,256,tmpName,params.jobid);
			hist2.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"monSumRHAC_%05d.csv");
			snprintf(csvName,256,tmpName,params.jobid);
			save1DRootCSV_float(csvName, hist2, nBins);
			
			snprintf(tmpName,256,dirName,"monSumSP_%05d.root");
			snprintf(rootName,256,tmpName,params.jobid);
			hist3.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"monSumSP_%05d.csv");
			snprintf(csvName,256,tmpName,params.jobid);
			save1DRootCSV_float(csvName, hist3, nBins);
			
			snprintf(tmpName,256,dirName,"monSumDS_%05d.root");
			snprintf(rootName,256,tmpName,params.jobid);
			hist4.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"monSumDS_%05d.csv");
			snprintf(csvName,256,tmpName,params.jobid);
			save1DRootCSV_float(csvName, hist4, nBins);
			
			snprintf(tmpName,256,dirName,"monSumRH_%05d.root");
			snprintf(rootName,256,tmpName,params.jobid);
			hist5.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"monSumRH_%05d.csv");
			snprintf(csvName,256,tmpName,params.jobid);
			save1DRootCSV_float(csvName, hist5, nBins);
			} break;
		case 11: {			
			// Save the root files
			snprintf(tmpName,256,dirName,"totalPMT1_%05d.root");
			snprintf(rootName,256,tmpName,params.jobid);
			hist1.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"totalPMT2_%05d.root");
			snprintf(rootName,256,tmpName,params.jobid);
			hist2.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"totalPMTC_%05d.root");
			snprintf(rootName,256,tmpName,params.jobid);
			hist3.SaveAs(rootName);
			// Save an output .csv file
			snprintf(tmpName,256,dirName,"totalByPMT_%05d.csv");
			snprintf(csvName,256,tmpName,params.jobid);
			save3DRootCSV(csvName,hist1,hist2,hist3,nBins);
		} break;
		case 12: {
			// Save the root files
			snprintf(tmpName,256,dirName,"coincPE1_S.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist1.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"coincPE2_S.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist2.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"coincPEC_S.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist3.SaveAs(rootName);
			// Save an output .csv file
			snprintf(tmpName,256,dirName,"coincPE_S.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			save3DRootCSV(csvName,hist1,hist2,hist3,nBins);
			// Save the root files
			snprintf(tmpName,256,dirName,"coincPE1_L.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist4.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"coincPE2_L.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist5.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"coincPEC_L.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			hist6.SaveAs(rootName);
			// Save an output .csv file
			snprintf(tmpName,256,dirName,"coincPE_L.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			save3DRootCSV(csvName,hist4,hist5,hist6,nBins);
		} break;
		default: {printf("Invalid setting %d for function 2!\n",setting); } break;
	}		
}

void functionality_background(std::vector<int> runList, params_in params, int setting) {
// For creating background .csv files. 	
// BE SUPER CAREFUL HERE -- There's no tagbit for heights in these.
// Have to pray we've properly separated beamOnBkg so hardcoding works.
// For future analyzers -- you'll probably have to change the heights.				

	const char* fillLoc  = std::getenv("FILL_LOC");  // MAD fill fit time const.
	const char* traceLoc = std::getenv("TRACE_LOC"); // Exp fill time const.
	const char* saveLoc  = std::getenv("SAVE_LOC");  // Save location for ROOT drawings
	
	char* fName;
	for (auto rIt = runList.begin(); rIt < runList.end(); rIt++) {
	
		int runNo = *rIt;
		fName = getFileName(runNo);
		std::ifstream infile(fName);
		if (infile.good()) {
			Run runMCS1(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_0");
			Run runMCS2(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_1");
			EMS runEMS(runNo, fName, "tems");
		
			if (runNo >= 9020) { runMCS2.setMCSOff(10);	}
			
			switch(setting) {
				case 1: {beamOffBkg(&runMCS1, &runMCS2, &runEMS); } break;
				case 2: {bkgRun1Bkg(&runMCS1, &runMCS2, &runEMS); } break;
				case 3: {bkgRun2Bkg(&runMCS1, &runMCS2, &runEMS); } break;
				case 4: {beamOnBkg(&runMCS1,&runMCS2, &runEMS);   } break;
				case 5: {productionBkg(&runMCS1,&runMCS2,&runEMS); } break;
				case 6: {beamOffMonitor(&runMCS1, &runMCS2, &runEMS); } break;
				default:{
					printf("Invalid setting %d for function 'Background'!\n",setting); break;
				} 
			}
		}
		fflush(stdout);
	}
					
					
}

void functionality_coinc(std::vector<int> runList, params_in params, int setting) {
// These functions are designed to study the timing information of pulse chains
// Like summing, it requires some ROOT histograms as the basic info.	
	const char* fillLoc  = std::getenv("FILL_LOC");  // MAD fill fit time const.
	const char* traceLoc = std::getenv("TRACE_LOC"); // Exp fill time const.
	const char* saveLoc  = std::getenv("SAVE_LOC");  // Save location for ROOT drawings
	
	// Cases 1 - 4 are timing + one 2D plot.
	int nTiming1D;
	int nTiming2D1;
	int nTiming2D2;
	switch (setting) {
		case 1: { // Cases 1-4 are various PMT timing and PHS 
			nTiming1D = 64000;
			nTiming2D1 = 75;
			nTiming2D2 = 75;
		} break;
		case 2: {
			nTiming1D = 64000;
			nTiming2D1 = 75;
			nTiming2D2 = 75;
		} break;
		case 3: {
			nTiming1D = 64000;
			nTiming2D1 = 75;
			nTiming2D2 = 75;
		} break;
		case 4: {
			nTiming1D = 64000;
			nTiming2D1 = 75;
			nTiming2D2 = 75;
		} break;
		case 5: {
			nTiming1D = 5000;
			nTiming2D1 = 5;
			nTiming2D2 = 5;
		} break;
		case 6: {
			nTiming1D = 5000;
			nTiming2D1 = 5;
			nTiming2D2 = 5;
		} break;
		case 7: {
			nTiming1D = 5000;
			nTiming2D1 = 1000;
			nTiming2D2 = 75;
		} break;
		case 8: {
			nTiming1D = 5000;
			nTiming2D1 = 1000;
			nTiming2D2 = 75;
		} break;
		case 9: {
			nTiming1D = 5000;
			nTiming2D1 = 5;
			nTiming2D2 = 5;
		} break;
		default : { 
			printf("Invalid setting %d for function 'Coinc' (4)!\n",setting);
		}
	}
	
	TH1D pTiming1("pTiming1", "pTiming1", nTiming1D+1,0,nTiming1D); 
	TH1D pTiming2("pTiming2", "pTiming2", nTiming1D+1,0,nTiming1D);
	TH1D pTiming3("pTiming3", "pTiming3", nTiming1D+1,0,nTiming1D);
	TH1D pTiming4("pTiming4", "pTiming4", nTiming1D+1,0,nTiming1D);
	
	// phsHits is an integer!
	TH2I phsHits("phsHits", "phsHits",nTiming2D1+1,0,nTiming2D1,nTiming2D2+1,0,nTiming2D2); // Pulse height spectrum
	// Pulse timing spectrum (2D -- time and number of hits), doubles.
	TH2D pT2D1("pTiming2D1","pTiming2D1",nTiming2D1+1,1,nTiming2D1,nTiming2D2+1,1,nTiming2D2);
	TH2D pT2D2("pTiming2D2","pTiming2D2",nTiming2D1+1,1,nTiming2D1,nTiming2D2+1,1,nTiming2D2);
	TH2D pT2D3("pTiming2D3","pTiming2D3",nTiming2D1+1,1,nTiming2D1,nTiming2D2+1,1,nTiming2D2);
	TH2D pT2D4("pTiming2D4","pTiming2D4",nTiming2D1+1,1,nTiming2D1,nTiming2D2+1,1,nTiming2D2);
		
	char* fName;
	for (auto rIt = runList.begin(); rIt < runList.end(); rIt++) {
	
		int runNo = *rIt;
		fName = getFileName(runNo);
		std::ifstream infile(fName);
		if (infile.good()) {
			Run runMCS1(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_0");
			Run runMCS2(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_1");
			EMS runEMS(runNo, fName, "tems");
		
			if (runNo >= 9020) { runMCS2.setMCSOff(10);	}
			// While cleaning this up, I found that I open this file?
			// FILE* conffile;
			// conffile = fopen("PMTRates.csv", "a");
			switch(setting) {
				// For note sake, swapping bool "background" to int "background" -- 0 = true, 1 = false, 2/3 new behavior
				case 1: { savePMTHits(&runMCS1, &runMCS2, &pTiming1,&pTiming2,&pTiming3,&pTiming4, &phsHits, 0);} break;
				case 2: { savePMTHits(&runMCS1, &runMCS2, &pTiming1,&pTiming2,&pTiming3,&pTiming4, &phsHits, 1);} break;
				case 3: { savePMTHits(&runMCS1, &runMCS2, &pTiming1,&pTiming2,&pTiming3,&pTiming4, &phsHits, 2);} break;
				case 4: { savePMTHits(&runMCS1, &runMCS2, &pTiming1,&pTiming2,&pTiming3,&pTiming4, &phsHits, 3);} break;
				case 5: { singlePMTHits(&runMCS1, &runMCS2, &pTiming1, &pTiming2); } break;
				case 6: { singlePMTHitsMoving(&runMCS1, &runMCS2, &pTiming1, &pTiming2); } break;
				case 7: { savePMTHits2D(&runMCS1, &runMCS2, &pT2D1,&pT2D2,&pT2D3,&pT2D4, true);} break;
				case 8: { savePMTHits2D(&runMCS1, &runMCS2, &pT2D1,&pT2D2,&pT2D3,&pT2D4, false);} break;
				case 9: { fillDetHitsMoving(&runMCS1, &runMCS2, &pTiming1);} break;
				default:{
					printf("Invalid setting %d for function 'Coinc' (4)!\n",setting);
				} 
			}
		}
		fflush(stdout);
	}
	
	// Initialize saving names
	char dirName[256];
	snprintf(dirName, 256, "%s/%%s", saveLoc);
	char rootName[256];
	char csvName[256];
	char tmpName[256];	
	// TODO: Clean this!
	switch(setting) {
		case 1: { 
			// Save histograms (ROOT)
			//std::string tmpName;
			snprintf(rootName,256,dirName,"pTiming1Bkg.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming1.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming2Bkg.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming2.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming3Bkg.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming3.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming4Bkg.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming4.SaveAs(tmpName);
			phsHits.Draw("COLZ");
			snprintf(rootName,256,dirName,"phsHitsBkg.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			phsHits.SaveAs(tmpName);
			
			int binNo = 0;
			
			// Convert histograms into .csv file
			
			snprintf(csvName,256,dirName,"phsHitsBkg.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* outfile = fopen(tmpName, "w");
			for (int ii = 1; ii <= nTiming2D1; ii++) {
				for (int jj = 1; jj <= nTiming2D2; jj++) {
					binNo = phsHits.GetBin(ii,jj);
					fprintf(outfile,"%d,%d:%d\n",(int)phsHits.GetBinContent(binNo),ii,jj);
				}
			}					
			fclose(outfile);
			// Save timing histos as csv
			snprintf(csvName,256,dirName,"pTiming1Bkg.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile1 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming2Bkg.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile2 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming3Bkg.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile3 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming4Bkg.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile4 = fopen(tmpName,"w");
			for (int ii = 1; ii <= nTiming1D; ii++) {
				fprintf(pOutfile1,"%d,%d\n",(int)pTiming1.GetBinContent(ii),ii);
				fprintf(pOutfile2,"%d,%d\n",(int)pTiming2.GetBinContent(ii),ii);
				fprintf(pOutfile3,"%d,%d\n",(int)pTiming3.GetBinContent(ii),ii);
				fprintf(pOutfile4,"%d,%d\n",(int)pTiming4.GetBinContent(ii),ii);
			}
			fclose(pOutfile1);
			fclose(pOutfile2);
			fclose(pOutfile3);
			fclose(pOutfile4);
			
			} break;
		case 2: {
			// Save histograms (ROOT)
			snprintf(rootName,256,dirName,"pTiming1.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming1.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming2.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming2.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming3.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming3.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming4.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming4.SaveAs(tmpName);
			phsHits.Draw("COLZ");
			snprintf(rootName,256,dirName,"phsHits.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			phsHits.SaveAs(tmpName);
			
			int binNo = 0;
			
			// Convert histograms into .csv file
			snprintf(csvName,256,dirName,"phsHits.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* outfile = fopen(tmpName, "w");
			for (int ii = 1; ii <= nTiming2D1; ii++) {
				for (int jj = 1; jj <= nTiming2D2; jj++) {
					binNo = phsHits.GetBin(ii,jj);
					fprintf(outfile,"%d,%d:%d\n",(int)phsHits.GetBinContent(binNo),ii,jj);
				}
			}					
			fclose(outfile);
			// Save timing histos as csv
			snprintf(csvName,256,dirName,"pTiming1.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile1 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming2.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile2 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming3.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile3 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming4.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile4 = fopen(tmpName,"w");
			for (int ii = 1; ii <= nTiming1D; ii++) {
				fprintf(pOutfile1,"%d,%d\n",(int)pTiming1.GetBinContent(ii),ii);
				fprintf(pOutfile2,"%d,%d\n",(int)pTiming2.GetBinContent(ii),ii);
				fprintf(pOutfile3,"%d,%d\n",(int)pTiming3.GetBinContent(ii),ii);
				fprintf(pOutfile4,"%d,%d\n",(int)pTiming4.GetBinContent(ii),ii);
			}
			fclose(pOutfile1);
			fclose(pOutfile2);
			fclose(pOutfile3);
			fclose(pOutfile4);
			} break;
		case 3: {
			// Save histograms (ROOT)
			snprintf(rootName,256,dirName,"pTiming1HoldS.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming1.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming2HoldS.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming2.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming3HoldS.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming3.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming4HoldS.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming4.SaveAs(tmpName);
			phsHits.Draw("COLZ");
			snprintf(rootName,256,dirName,"phsHitsHoldS.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			phsHits.SaveAs(tmpName);
			
			int binNo = 0;
			
			// Convert histograms into .csv file
			snprintf(csvName,256,dirName,"phsHitsHoldS.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* outfile = fopen(tmpName, "w");
			for (int ii = 1; ii <= nTiming2D1; ii++) {
				for (int jj = 1; jj <= nTiming2D2; jj++) {
					binNo = phsHits.GetBin(ii,jj);
					fprintf(outfile,"%d,%d:%d\n",(int)phsHits.GetBinContent(binNo),ii,jj);
				}
			}					
			fclose(outfile);
			// Save timing histos as csv
			snprintf(csvName,256,dirName,"pTiming1HoldS.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile1 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming2HoldS.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile2 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming3HoldS.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile3 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming4HoldS.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile4 = fopen(tmpName,"w");
			for (int ii = 1; ii <= nTiming1D; ii++) {
				fprintf(pOutfile1,"%d,%d\n",(int)pTiming1.GetBinContent(ii),ii);
				fprintf(pOutfile2,"%d,%d\n",(int)pTiming2.GetBinContent(ii),ii);
				fprintf(pOutfile3,"%d,%d\n",(int)pTiming3.GetBinContent(ii),ii);
				fprintf(pOutfile4,"%d,%d\n",(int)pTiming4.GetBinContent(ii),ii);
			}
			fclose(pOutfile1);
			fclose(pOutfile2);
			fclose(pOutfile3);
			fclose(pOutfile4);
			} break;
		case 4 : { 
		// Save histograms (ROOT)
			snprintf(rootName,256,dirName,"pTiming1Unl.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming1.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming2Unl.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming2.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming3Unl.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming3.SaveAs(tmpName);
			snprintf(rootName,256,dirName,"pTiming4Unl.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			pTiming4.SaveAs(tmpName);
			phsHits.Draw("COLZ");
			snprintf(rootName,256,dirName,"phsHitsUnl.root%d");
			snprintf(tmpName,256,rootName,params.jobid);
			phsHits.SaveAs(tmpName);
			
			int binNo = 0;
			
			// Convert histograms into .csv file
			snprintf(csvName,256,dirName,"phsHitsUnl.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* outfile = fopen(tmpName, "w");
			for (int ii = 1; ii <= nTiming2D1; ii++) {
				for (int jj = 1; jj <= nTiming2D2; jj++) {
					binNo = phsHits.GetBin(ii,jj);
					fprintf(outfile,"%d,%d:%d\n",(int)phsHits.GetBinContent(binNo),ii,jj);
				}
			}					
			fclose(outfile);
			// Save timing histos as csv
			snprintf(csvName,256,dirName,"pTiming1Unl.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile1 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming2Unl.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile2 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming3Unl.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile3 = fopen(tmpName,"w");
			snprintf(csvName,256,dirName,"pTiming4Unl.csv%d");
			snprintf(tmpName,256,csvName,params.jobid);
			FILE* pOutfile4 = fopen(tmpName,"w");
			for (int ii = 1; ii <= nTiming1D; ii++) {
				fprintf(pOutfile1,"%d,%d\n",(int)pTiming1.GetBinContent(ii),ii);
				fprintf(pOutfile2,"%d,%d\n",(int)pTiming2.GetBinContent(ii),ii);
				fprintf(pOutfile3,"%d,%d\n",(int)pTiming3.GetBinContent(ii),ii);
				fprintf(pOutfile4,"%d,%d\n",(int)pTiming4.GetBinContent(ii),ii);
			}
			fclose(pOutfile1);
			fclose(pOutfile2);
			fclose(pOutfile3);
			fclose(pOutfile4);
			} break;
		case 5: {
			snprintf(tmpName,256,dirName,"pmt1DT.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pTiming1.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"pmt2DT.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pTiming2.SaveAs(rootName);		
			} break;
		case 6: {
			snprintf(tmpName,256,dirName,"pmt1DT.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pTiming1.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"pmt2DT.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pTiming2.SaveAs(rootName);					
		} break;
		case 7: {// Save histograms (ROOT)
			snprintf(tmpName,256,dirName,"pT2D1_bkg.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pT2D1.Draw("COLZ");
			pT2D1.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"pT2D2_bkg.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pT2D2.Draw("COLZ");
			pT2D2.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"pT2D3_bkg.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pT2D3.Draw("COLZ");
			pT2D3.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"pT2D4_bkg.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pT2D4.Draw("COLZ");
			pT2D4.SaveAs(rootName);
			
			int binNo = 0;
			
			// Save timing histos as csv
			snprintf(tmpName,256,dirName,"pT2D1_bkg.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			FILE* pOut2DB1 = fopen(csvName,"w");
			snprintf(tmpName,256,dirName,"pT2D2_bkg.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			FILE* pOut2DB2 = fopen(csvName,"w");
			snprintf(tmpName,256,dirName,"pT2D3_bkg.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			FILE* pOut2DB3 = fopen(csvName,"w");
			snprintf(tmpName,256,dirName,"pT2D4_bkg.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			FILE* pOut2DB4 = fopen(csvName,"w");
			for (int ii = 1; ii <= nTiming2D1; ii++) {
				for (int jj = 1; jj <= nTiming2D2; jj++) {
					binNo = pT2D1.GetBin(ii,jj);
					fprintf(pOut2DB1,"%d,%d,%d\n",(int)pT2D1.GetBinContent(binNo),ii,jj);
					fprintf(pOut2DB2,"%d,%d,%d\n",(int)pT2D2.GetBinContent(binNo),ii,jj);
					fprintf(pOut2DB3,"%d,%d,%d\n",(int)pT2D3.GetBinContent(binNo),ii,jj);
					fprintf(pOut2DB4,"%d,%d,%d\n",(int)pT2D4.GetBinContent(binNo),ii,jj);
				}
			}
			fclose(pOut2DB1);
			fclose(pOut2DB2);
			fclose(pOut2DB3);
			fclose(pOut2DB4);
			} break;
		case 8: {
			snprintf(tmpName,256,dirName,"pT2D1.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pT2D1.Draw("COLZ");
			pT2D1.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"pT2D2.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pT2D2.Draw("COLZ");
			pT2D2.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"pT2D3.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pT2D3.Draw("COLZ");
			pT2D3.SaveAs(rootName);
			snprintf(tmpName,256,dirName,"pT2D4.root%d");
			snprintf(rootName,256,tmpName,params.jobid);
			pT2D4.Draw("COLZ");
			pT2D4.SaveAs(rootName);
			
			int binNo = 0;
			// Save timing histos as csv
			snprintf(tmpName,256,dirName,"pT2D1.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			FILE* pOut2D1 = fopen(csvName,"w");
			snprintf(tmpName,256,dirName,"pT2D2.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			FILE* pOut2D2 = fopen(csvName,"w");
			snprintf(tmpName,256,dirName,"pT2D3.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			FILE* pOut2D3 = fopen(csvName,"w");
			snprintf(tmpName,256,dirName,"pT2D4.csv%d");
			snprintf(csvName,256,tmpName,params.jobid);
			FILE* pOut2D4 = fopen(csvName,"w");
			for (int ii = 1; ii <= nTiming2D1; ii++) {
				for (int jj = 1; jj <= nTiming2D2; jj++) {
					binNo = pT2D1.GetBin(ii,jj);
					fprintf(pOut2D1,"%d,%d,%d\n",(int)pT2D1.GetBinContent(binNo),ii,jj);
					fprintf(pOut2D2,"%d,%d,%d\n",(int)pT2D2.GetBinContent(binNo),ii,jj);
					fprintf(pOut2D3,"%d,%d,%d\n",(int)pT2D3.GetBinContent(binNo),ii,jj);
					fprintf(pOut2D4,"%d,%d,%d\n",(int)pT2D4.GetBinContent(binNo),ii,jj);
				}
			}
			fclose(pOut2D1);
			fclose(pOut2D2);
			fclose(pOut2D3);
			fclose(pOut2D4);
			} break;
		case 9: {
			snprintf(csvName,256,"fillDT_%d.root",params.coincMode);
			snprintf(rootName,256, dirName,csvName);
			pTiming1.SaveAs(rootName);
		} break;
		default : {printf("Invalid setting %d for function 'Coinc' (4)!\n",setting); } break;
	}

}

void functionality_like(std::vector<int> runList, params_in params, int setting) {
// These were some tests I did with likelihood analyzer tests. 
// Someone smarter than I am could probably do stuff with them, I haven't 
// messed around with them enough to get a good conclusion.
//TODO: MAKE GO (I've commented it out to get compilation)
	// Function Extraction
	const char* fillLoc  = std::getenv("FILL_LOC");  // MAD fill fit time const.
	const char* traceLoc = std::getenv("TRACE_LOC"); // Exp fill time const.
	const char* saveLoc  = std::getenv("SAVE_LOC");  // Save location for ROOT drawings
	
	//// Likelihood study histograms
	//TH2D lTMax("lTMax","lTMax",1001,0,10000,100,0,500);
	//TH2D lNPE("lNPE", "lNPE",50,0,49, 100,0,500);
	
	char* fName;
	for (auto rIt = runList.begin(); rIt < runList.end(); rIt++) {
	
		int runNo = *rIt;
		fName = getFileName(runNo);
		std::ifstream infile(fName);
		if (infile.good()) {
			Run runMCS1(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_0");
			Run runMCS2(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_1");
			EMS runEMS(runNo, fName, "tems");
		
			if (runNo >= 9020) { runMCS2.setMCSOff(10);	}
			
			switch(setting) {
				//case 1: { likelihoodStudy(&runMCS1, &runMCS2, &lNPE, &lTMax); } break;
				//case 2: { multiCoincCheck(&runMCS1,&runMCS2,&numPh,&coinLen,&phByLen);} break;
				default:{
					printf("Invalid setting %d for function 'Likelihood'!\n",setting);
				} 
			}
		}
		fflush(stdout);
	}
}

void functionality_fastPMT(std::vector<int> runList, params_in params, int setting) {
// These were designed to look specifically at "fast" pmt coincidence events
// The different cases here are different types of runs/sections of running.

	// Function Extraction
	const char* fillLoc  = std::getenv("FILL_LOC");  // MAD fill fit time const.
	const char* traceLoc = std::getenv("TRACE_LOC"); // Exp fill time const.
	const char* saveLoc  = std::getenv("SAVE_LOC");  // Save location for ROOT drawings
	
	int bin2d = 101; // bin divider
	int photons1 = 100; // photon number
	int photons2 = 100; // I think I had some variable in the asym number?
	int time = 10000; // number of timesteps
	
	TH2D time2d("photons","photons",bin2d,0,photons1,bin2d,0,time);
	TH2D asym2d("asByL","asByL",bin2d,0,photons1,bin2d,0,photons2);
	
	char* fName;
	for (auto rIt = runList.begin(); rIt < runList.end(); rIt++) {
	
		int runNo = *rIt;
		fName = getFileName(runNo);
		std::ifstream infile(fName);
		if (infile.good()) {
			Run runMCS1(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_0");
			Run runMCS2(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_1");
			EMS runEMS(runNo, fName, "tems");
		
			if (runNo >= 9020) { runMCS2.setMCSOff(10);	}
			
			switch(setting) {
				case 1: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 0);} break;
				case 2: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 1);} break;
				case 3: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 2);} break;
				case 4: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 3);} break;
				case 5: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 4);} break;
				case 6: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 5);} break;
				case 7: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 6);} break;
				case 8: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 7);} break;
				default:{
					printf("Invalid setting %d for function 'Fast Coinc'!\n",setting);
				} break;
			} 
		}
		fflush(stdout);
	}
	// Now initialize saving names:
	// We do the same thing for all of these, so it's just a name switch.
	char dirName[256];
	snprintf(dirName,256, "%s/%%s", saveLoc);
	char tmpNameT[256];
	char tmpNameA[256];
	char rootNameT[256];
	char rootNameA[256];
	char csvNameT[256];
	char csvNameA[256];	
	
	switch (setting) { 
		case 1: {
			snprintf(rootNameT,256,dirName,"hitsByLenBkg38.root%d");
			snprintf(rootNameA,256,dirName,"asymByLenBkg38.root%d");
			snprintf(csvNameT,256,dirName,"hitsByLenBkg38.csv%d");
			snprintf(csvNameA,256,dirName,"asymByLenBkg38.csv%d");
		} break;
		case 2: {
			snprintf(rootNameT,256,dirName,"hitsByLenBkg49.root%d");
			snprintf(rootNameA,256,dirName,"asymByLenBkg49.root%d");
			snprintf(csvNameT,256,dirName,"hitsByLenBkg49.csv%d");
			snprintf(csvNameA,256,dirName,"asymByLenBkg49.csv%d");
		} break;
		case 3: {
			snprintf(rootNameT,256,dirName,"hitsByLenBkg25.root%d");
			snprintf(rootNameA,256,dirName,"asymByLenBkg25.root%d");
			snprintf(csvNameT,256,dirName,"hitsByLenBkg25.csv%d");
			snprintf(csvNameA,256,dirName,"asymByLenBkg25.csv%d");
		} break;
		case 4: {
			snprintf(rootNameT,256,dirName,"hitsByLenBkg1.root%d");
			snprintf(rootNameA,256,dirName,"asymByLenBkg1.root%d");
			snprintf(csvNameT,256,dirName,"hitsByLenBkg1.csv%d");
			snprintf(csvNameA,256,dirName,"asymByLenBkg1.csv%d");
		} break;
		case 5: {
			snprintf(rootNameT,256,dirName,"hitsByLenBkg.root%d");
			snprintf(rootNameA,256,dirName,"asymByLenBkg.root%d");
			snprintf(csvNameT,256,dirName,"hitsByLenBkg.csv%d");
			snprintf(csvNameA,256,dirName,"asymByLenBkg.csv%d");
		} break;
		case 6: {
			snprintf(rootNameT,256,dirName,"hitsByLen.root%d");
			snprintf(rootNameA,256,dirName,"asymByLen.root%d");
			snprintf(csvNameT,256,dirName,"hitsByLen.csv%d");
			snprintf(csvNameA,256,dirName,"asymByLen.csv%d");
		} break;
		case 7: {
			snprintf(rootNameT,256,dirName,"hitsByLenHold.root%d");
			snprintf(rootNameA,256,dirName,"asymByLenHold.root%d");
			snprintf(csvNameT,256,dirName,"hitsByLenHold.csv%d");
			snprintf(csvNameA,256,dirName,"asymByLenHold.csv%d");
		} break;
		case 8: {
			snprintf(rootNameT,256,dirName,"hitsByLenUnl.root%d");
			snprintf(rootNameA,256,dirName,"asymByLenUnl.root%d");
			snprintf(csvNameT,256,dirName,"hitsByLenUnl.csv%d");
			snprintf(csvNameA,256,dirName,"asymByLenUnl.csv%d");
		} break;
		default:{
			printf("Invalid setting %d for function 'Fast Coinc'!\n",setting);
		} 
	}
	snprintf(tmpNameT,256,rootNameT,params.jobid);
	time2d.SaveAs(tmpNameT);
	snprintf(tmpNameA,256,rootNameA,params.jobid);
	asym2d.SaveAs(tmpNameA);
	
	// Hardcode save this as csv
	// TODO functionalize this
	int binNo = 0;
	// Convert histograms into .csv file
	snprintf(tmpNameT,256,csvNameT,params.jobid);
	FILE* outfileH = fopen(tmpNameT, "w");
	for (int ii = 1; ii <= bin2d; ii++) {
		for (int jj = 1; jj <= bin2d; jj++) {
			binNo = time2d.GetBin(ii,jj);
			fprintf(outfileH,"%d,%d:%d\n",(int)time2d.GetBinContent(binNo),ii,jj);
		}
	}					
	fclose(outfileH);
	binNo = 0;
	// Convert histograms into .csv file, same for asym
	snprintf(tmpNameA,256,csvNameA,params.jobid);
	FILE* outfileA = fopen(tmpNameA, "w");
	for (int ii = 1; ii <= bin2d; ii++) {
		for (int jj = 1; jj <= bin2d; jj++) {
			binNo = asym2d.GetBin(ii,jj);
			fprintf(outfileA,"%d,%d:%d\n",(int)asym2d.GetBinContent(binNo),ii,jj);
		}
	}					
	fclose(outfileA);	
}

int main(int argc, const char** argv) {

	char fName[256];

	using namespace std::placeholders;
	if((argc > 1 && !strcmp(argv[1], "help")) || argc < 8) {
		printf("\nUsage: ./Analyzer 'QUERY_Runs' coincWindow peSumWindow peSum coincMode func set (ID)\n");
		
		printf("\nQUERY_Runs should be written as:\n");
		printf("   A file ");
		printf("   (Or just a run number if you want...)\n");
		
		printf("\nDefault Values (From 2017-2018): \n");
		printf("   coincWindow = 50\n");
		printf("   peSumWindow = 500\n");
		printf("   peSum = 8\n");
		printf("\nPossible Coincidence Modes: \n");
		printf("   1: Fixed Window (Low threshold)\n");
		printf("   2: Fixed Window (High threshold)\n");
		printf("   3: Moving Window (Low threshold)\n");
		printf("   4: Moving Window (High threshold)\n");
		printf("   5: Arbitrary PE (Low threshold)\n");
		printf("   6: Arbitrary PE (High Threshold)\n");
		printf("   func and set tables below. \n");
				
		printf("\nFunctionality: \n");
		printf("  -1 --> Testing new functions\n");
		printf("   1 --> calculate normalized values by dip\n");
		printf("         Setting: 1 --> use normByNDip\n");
		printf("         Setting: 2 --> use normByNDipSing\n");
		printf("         Setting: 3 --> use normByNDipCoinc\n");
		printf("   2 --> sum histograms\n");
		printf("         Setting: 1 --> use dipSummerS\n");
		printf("         Setting: 2 --> use dipSummerL\n");
		printf("         Setting: 3 --> use fillSummer\n");
		printf("         Setting: 4 --> use fitLastDip\n");
		printf("         Setting: 5 --> use fitFill\n");
		printf("         Setting: 6 --> use acSummerCoinc\n");
		printf("         Setting: 7 --> use longHoldSummer\n");
		printf("         Setting: 8 --> use bkgSummer\n");
		printf("         Setting: 9 --> use peak1SummerS\n");
		printf("   3 --> calculate background values\n");
		printf("         Setting: 1 --> use beamOffBkg\n");
		printf("         Setting: 2 --> use bkgRun1Bkg\n");
		printf("         Setting: 3 --> use bkgRun2Bkg\n");
		printf("         Setting: 4 --> use beamOnBkg\n");
		printf("         Setting: 5 --> use productionBkg\n");
		printf("   4 --> PMT and Coincidence analysis\n");
		printf("         Setting: 1 --> use findPhotons (background) \n");
		printf("         Setting: 2 --> use findPhotons (foreground) \n");
		printf("         Setting: 3 --> use nextPhoton \n");
		printf("         Setting: 4 --> use nextPhotonMultiple \n");
		printf("         Setting: 5 --> use savePMTHits2D (background) \n");
		printf("         Setting: 6 --> use savePMTHits2D (foreground) \n");
		printf("         Setting: 7 --> use f");
		printf("   5 --> Likelihood filtering studies\n");
		printf("         Setting: 1 --> use likelihoodStudy \n");
		printf("         Setting: 2 --> use multiCoincCheck \n");
		printf("   Other func and set tables not yet incorporated.\n");
				
		return 1;
	}
	
	// Trying to get 2017 and 2018 to both work.
	char* fileLoc17 = std::getenv("RUNDIR_2017");
	char* fileLoc18 = std::getenv("RUNDIR_2018");
	char* fileLoc19 = std::getenv("RUNDIR_2019");
	
	const char* fillLoc  = std::getenv("FILL_LOC");  // MAD fill fit time const.
	const char* traceLoc = std::getenv("TRACE_LOC"); // Exp fill time const.
	const char* saveLoc  = std::getenv("SAVE_LOC");  // Save location for ROOT drawings
	
	if ((fileLoc17) || (fileLoc18) || (fileLoc19)) {
		if (fileLoc17) {
			printf("Using 2017 run location %s\n", fileLoc17);
		}
		if (fileLoc18) {
			printf("Using 2018 run location %s\n", fileLoc18);
		}
		if (fileLoc19) {
			printf("Using 2019 run location %s\n", fileLoc19);
		}
	} else {
		printf("ERROR: Environmental variable RUNDIR not set! \n");
		return 0;
	}
	 
	// Initialize Variables
	std::string query = argv[1];
	printf("This was the query sent:\n%s\n", query.c_str());
	
	params_in params; // Coincidence Settings -- Set up to go to list
	// TODO: Diversify, to make this a loaded file (so that we can do scans)
	params.coincWindow = atoi(argv[2]);
	params.peSumWindow = atoi(argv[3]);
	params.peSum       = atoi(argv[4]);
	params.coincMode   = atoi(argv[5]);
	
	// int coincWindow   = atoi(argv[2]);
	// int peSumWindow   = atoi(argv[3]);
	//int peSum         = atoi(argv[4]);
	//int coincMode     = atoi(argv[5]);
	int functionality = atoi(argv[6]);
	int setting       = atoi(argv[7]);
	
	const char* jobidTmp = argv[8]; // if we're gnuparallel-ing things
	//int jobid;
	if (jobidTmp != NULL) { // Check -- if we didn't parallelize, then set ID as 0
		params.jobid = atoi(argv[8]);
		query += std::to_string(params.jobid);
	} else {
		params.jobid = 0;
	}
	
	// TODO: maybe a lookup table (combine with the intro check)
	if ((params.coincMode > 999) && (functionality !=4)){
		printf("Invalid coincidence mode specified! \n");
		return 1;
	}
	if (functionality > 15 || setting > 15) {
		printf("Invalid functionality or setting! Check your table here...\n");
		return 1;
	}
	printf("Using functionality %d,%d\n",functionality,setting);
	//------------------------------------------------------------------
	// Create ROOT histograms (and other initializations) 
	// for our different functions!
	//
	// The ones we actually use vary based on which run we're doing.
	// However, we can't declare them in a loop if we want to run
	// across everything, so they're stuck being declared here.
	//
	// Like all of these things should be removed
	//------------------------------------------------------------------
		
	//double endFillTime;
	//std::vector<int> numHGXhits;
	//std::vector<double> beamHits;
	//std::vector<double> sCts;
	//std::vector<double> lCts;
	
	
			
	//// Sums of singles dips 
	//TH1D singleBkg("singleBkg", "singleBkg", 1000,-500,500);

	//TH1D dipSumACS("summedDipACS","summedDipAC",3100,0,310);
	//TH1D dipSumACL("summedDipACL","summedDipAC",3100,0,310);
	
	//// Sums of all data (300 + 50 + 20 + 210 + 50 = 630)
	//// And for a long hold (300 + 50 + 1550 + 210 + 50 = 2160
	//TH1D totalPMT1("summedPMT1","summedPMT1",2160,0,2160);
	//TH1D totalPMT2("summedPMT2","summedPMT2",2160,0,2160);
	//TH1D totalCoin("summedPMTC","summedPMTC",2160,0,2160);
	
	//TH1D dayBkg("summedBkg","summedBkg",10000,0,1000);
	
	//// Background sums
	//TH1D bkgSum1("summedBkg1","summedBkg1", 51,0,50);
	//TH1D bkgSum2("summedBkg2","summedBkg2", 51,0,50);
	//// PMT hits information
	//int nTB = 50000; // number of timing bins
	//int phsTB = 75;  // number of pulse height bins
	
	//// Pulse timing spectrum for each combination of PMT and initial PMT hit
	//TH1D pTiming1("pTiming1", "pTiming1", nTB+1,0,nTB); 
	//TH1D pTiming2("pTiming2", "pTiming2", nTB+1,0,nTB);
	//TH1D pTiming3("pTiming3", "pTiming3", nTB+1,0,nTB);
	//TH1D pTiming4("pTiming4", "pTiming4", nTB+1,0,nTB);
	//TH2I phsHits("phsHits", "phsHits", phsTB,1,phsTB,phsTB,1,phsTB); // Pulse height spectrum
	
	//// Pulse timing spectrum for each PMT's spacing
	//char pmtName [16];
	//if ((functionality == 4) && (setting ==7)) {
		//sprintf(pmtName, "pmt_%d_DT", coincMode);
	//} else {
		//sprintf(pmtName, "pmt1DT");
	//}
	//TH1D pmt1DT(pmtName, pmtName, 5001,0,5000);
	//TH1D pmt2DT("pmt2DT", "pmt2DT", 1001,0,1000);
	
	//// Pulse timing spectrum (2D -- time and number of hits)
	//TH2D pT2D1("pTiming2D1","pTiming2D1",1001,0,1000,phsTB,1,phsTB);
	//TH2D pT2D2("pTiming2D2","pTiming2D2",1001,0,1000,phsTB,1,phsTB);
	//TH2D pT2D3("pTiming2D3","pTiming2D3",1001,0,1000,phsTB,1,phsTB);
	//TH2D pT2D4("pTiming2D4","pTiming2D4",1001,0,1000,phsTB,1,phsTB);
	
	//// Likelihood study histograms
	//TH2D lTMax("lTMax","lTMax",1001,0,10000,100,0,500);
	//TH2D lNPE("lNPE", "lNPE",50,0,49, 100,0,500);
	
	////TH2D TH1D* nPh, TH1D* coinL, TH2D* phByL
	//TH2D numPh("numPh", "numPh", 2601,0,2600,76,0,75);
	//TH2D coinLen("coinL","coinL",2601,0,2600,401,0,400);
	//TH2D phByLen("phByL","phByL",76,0,75,4001,0,4000);
	//// want 150 bins of 50 ns each so 7500
	//printf("Loaded ROOT hists\n");
	//------------------------------------------------------------------
	// Load run data from .root file. 
	// Remember that we need to set the global RUNDIR variable
	//------------------------------------------------------------------
	std::istringstream iss(query); // convert argv[1] into stringstream
	std::string token;
	
	//int nRuns = 0; // Run counter
	/*std::vector<int> rList;
	while(std::getline(iss, token, ',')) { // This should load a file
		int runNo = atoi(token.c_str());
		if (runNo > 0) { // If we just entered a run or list of runs
			rList.push_back(runNo);
		} else {
			std::ifstream run_file(query.c_str()); // Treats the query as a file
			if (run_file.is_open()) { // Check if the run_file works
				printf("Able to load file:\n");
				std::string line;
				while(std::getline(run_file,line,',')) {
					rList.push_back(std::stoi(line));
					printf("%05d,",std::stoi(line));
				}
			}
		}
	}*/
	
	
	printf("Loading run list\n");
	
	std::vector<int> rList = load_run_list(query);
	if (rList.size() == 0) { 
		printf("Error! Unable to load any runs!\n");
	}
	// Special Case: We're NOT parallelizing and running one run
	if ((jobidTmp == NULL) && (rList.size() == 1)){
		params.jobid = rList.at(0); // Just load the one run
	}
		
	switch (functionality) {
		case -1: { // Debug Panel! Not a separate function.
			char* fName;
			for (auto rIt = rList.begin(); rIt < rList.end(); rIt++) {
				int runNo = *rIt;
				fName = getFileName(runNo);
				std::ifstream infile(fName);
				if (infile.good()) {
					Run runMCS1(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_0");
					Run runMCS2(params.coincWindow, params.peSumWindow, params.peSum, runNo, params.coincMode, fName, "tmcs_1");
					EMS runEMS(runNo, fName, "tems");
					// Put Debug Functions here.
					int runType = tagbitTagger(&runMCS1, &runMCS2, true); 
					int value = coincidenceCompare(&runMCS1, &runMCS2);
				}
			}
		} break;
		case 1: {
			functionality_extract(rList,params,setting);
		} break;
		case 2: {
			functionality_summing(rList,params,setting);
		} break;
		case 3: {
			functionality_background(rList,params,setting);
		} break;
		case 4: {
			functionality_coinc(rList,params,setting);
		} break;
		case 5: {
			functionality_like(rList,params,setting);
		} break;
		case 6: {
			functionality_fastPMT(rList,params,setting);
		} break;
	}
	
	return 0;
	
	//for (auto rIt = rList.begin(); rIt < rList.end(); rIt++) {
	////while(std::getline(iss, token, ',')) { // Trying to get rid of the dumb cat thing
		////int runNo = atoi(token.c_str());
		//int runNo = *rIt;
		//// Hard coding in the 2017/2018 year split at 9519
		//// Also hard coding in run 3412 as the beginning of 2017.
		//// 9519-9538 are in both directories BUT only in 2018 runlogdump.
		//if ((3412 <= runNo) && (runNo < 9519)) {
			//snprintf(fName, 256, "%s/processed_output_%%05d.root", fileLoc17);
		//} else if ((9519 <= runNo) && (runNo < 14732)) {
			//snprintf(fName, 256, "%s/processed_output_%%05d.root", fileLoc18);
		//} else if((14732 <= runNo) && (runNo < 17999)) {
			//snprintf(fName, 256, "%s/processed_output_%%05d.root", fileLoc19);
		//} else {
			//printf("Warning! Run %05d is not from calendar year 2017 or 2018 !\n", runNo);
			//continue;
		//}
				
		//std::string runBody = fName;
		//sprintf(fName, runBody.c_str(), runNo);
		//std::ifstream infile(fName);
		
		//if (infile.good()) {
			//Run runMCS1(coincWindow, peSumWindow, peSum, runNo, coincMode, fName, "tmcs_0");
			//Run runMCS2(coincWindow, peSumWindow, peSum, runNo, coincMode, fName, "tmcs_1");
			//EMS runEMS(runNo, fName, "tems");
			
			//// MCS2 changes by 10 at run 9020:
			//if (runNo >= 9020) {
				//runMCS2.setMCSOff(10);
			//}
					
			//// This is where we run our neccesary function. 
			//// It's a nested switch statement, sorry not sorry.	
			//switch (functionality) {
				//// Case -1 is the test function section.
				//case -1: {
					////std::vector<double> blah = getLXOn(&runMCS1, &runMCS2, &runEMS); 
					//int runType = tagbitTagger(&runMCS1, &runMCS2, true); 
					//int value = coincidenceCompare(&runMCS1, &runMCS2);
				//} break;
				
				//// Normalized run data
				//case 1 : {
					//switch(setting) {
						//case 1: {normNByDip(&runMCS1, &runMCS2, &runEMS, traceLoc);     } break;
						//case 2: {normNByDipSing(&runMCS1, &runMCS2, &runEMS, traceLoc); } break;
						//case 3: {normNByDipCoinc(&runMCS1, &runMCS2, &runEMS, traceLoc);} break;
						//case 4: {extractObservables(&runMCS1, &runMCS2, traceLoc);} break;
						//case 5: {activeCleaner(&runMCS1,&runMCS2,&runEMS,traceLoc);} break;
						//case 6: {fitBkgUnload(&runMCS1,&runMCS2,&runEMS);} break;
						//case 7: {normNByDipDet(&runMCS1,&runMCS2,&runEMS,traceLoc);} break;
						//default:{
							//printf("Invalid setting %d for function %d!\n",setting, functionality);
						//} 
					//} break;
				//// Histogram summing
				//} case 2 : { 
					///*switch(setting) { 
						//case 1: {dipSummerS(&runMCS1, &runMCS2, &dipSumS, 20.0); } break;
						//case 2: {dipSummerL(&runMCS1, &runMCS2, &dipSumL); } break;
						//case 3: {
							//for (int ii = 1; ii < 11; ii++) {
								//fitFillFree(&runMCS1,&runMCS2,ii);	
							//}} break;
							
						////case 4: {fitLastDip(&runMCS1, &runMCS2, &dipSumS);} break;
						//case 4: {writeFastEvts(&runMCS1,&runMCS2,&dipSumS);
								 //writeFastBkgs(&runMCS1,&runMCS2,&dayBkg);} break;
						////case 5: {fitFill(&runMCS1, &runMCS2, &fillSumBa);} break;
						//case 5: {
							//for (int ii = 1; ii < 11; ii++) {
								//fitFill(&runMCS1,&runMCS2,fillLoc,ii);	
							//}} break;
						//case 6: {acSummerCoinc(&runMCS1, &runMCS2, &acSumAll, 20.0);} break;
						//case 7: {longHoldSummer(&runMCS1, &runMCS2, &longPMT1, &longPMT2,&longCoinc);} break;
						//case 8: {bkgSummer(&runMCS1, &runMCS2, &bkgSum1,&bkgSum2);} break;
						//case 9: {peak1SummerS(&runMCS1, &runMCS2,&dipSumS,&dipSumACS);
								 //peak1SummerL(&runMCS1, &runMCS2,&dipSumL,&dipSumACL);} break;
						//case 10: {fillSummerD1(&runMCS1,350,&fillSumD1);
								  //fillSummerD2(&runMCS1,350,&fillSumD2);
								  //fillSummerOl(&runMCS1,350,&fillSumOl);
								  //fillSummerBa(&runMCS1,350,&fillSumBa);
								  //fillSummerSp(&runMCS2,350,&fillSumSp); } break;
						//case 11: {totalSummer(&runMCS1, &runMCS2, &totalPMT1, &totalPMT2, &totalCoin, 20.0);}break;
						//default: {
							//printf("Invalid setting %d for function %d!\n",setting, functionality);
						//} 
					//} break;
					//*/
				//// Background counting
				//} case 3 : {
					//switch(setting) { 
						//// BE SUPER CAREFUL HERE -- There's no tagbit for heights in these.
						//// Have to pray we've properly separated beamOnBkg so hardcoding works.
						//case 1: {beamOffBkg(&runMCS1, &runMCS2, &runEMS); } break;
						//case 2: {bkgRun1Bkg(&runMCS1, &runMCS2, &runEMS); } break;
						//case 3: {bkgRun2Bkg(&runMCS1, &runMCS2, &runEMS); } break;
						//case 4: {beamOnBkg(&runMCS1,&runMCS2, &runEMS);   } break;
						//case 5: {productionBkg(&runMCS1,&runMCS2,&runEMS); } break;
						//case 6: {beamOffMonitor(&runMCS1, &runMCS2, &runEMS); } break;
						//default:{
							//printf("Invalid setting %d for function %d!\n",setting, functionality);
						//} 
					//} break;
					
				//} case 4 : { 
					//FILE* conffile;
					//conffile = fopen("PMTRates.csv", "a");
					
					//switch(setting) {
						//// For note sake, swapping bool "background" to int "background" -- 0 = true, 1 = false, 2/3 new behavior
						//case 1: { savePMTHits(&runMCS1, &runMCS2, &pTiming1,&pTiming2,&pTiming3,&pTiming4, &phsHits, 0);} break;
						//case 2: { savePMTHits(&runMCS1, &runMCS2, &pTiming1,&pTiming2,&pTiming3,&pTiming4, &phsHits, 1);} break;
						//case 3: { savePMTHits(&runMCS1, &runMCS2, &pTiming1,&pTiming2,&pTiming3,&pTiming4, &phsHits, 2);} break;
						//case 4: { savePMTHits(&runMCS1, &runMCS2, &pTiming1,&pTiming2,&pTiming3,&pTiming4, &phsHits, 3);} break;
						//case 5: { singlePMTHits(&runMCS1, &runMCS2, &pmt1DT, &pmt2DT); } break;
						//case 6: { singlePMTHitsMoving(&runMCS1, &runMCS2, &pmt1DT, &pmt2DT); } break;
						//case 7: { savePMTHits2D(&runMCS1, &runMCS2, &pT2D1,&pT2D2,&pT2D3,&pT2D4, true);} break;
						//case 8: { savePMTHits2D(&runMCS1, &runMCS2, &pT2D1,&pT2D2,&pT2D3,&pT2D4, false);} break;
						//case 9: { fillDetHitsMoving(&runMCS1, &runMCS2, &pmt1DT);} break;
						//default:{
							//printf("Invalid setting %d for function %d!\n",setting, functionality);
						//} 
					//} break;
				//} case 5 : {
					//switch(setting) {
						//case 1: { likelihoodStudy(&runMCS1, &runMCS2, &lNPE, &lTMax); } break;
						//case 2: { multiCoincCheck(&runMCS1,&runMCS2,&numPh,&coinLen,&phByLen);} break;
						//default:{
							//printf("Invalid setting %d for function %d!\n",setting, functionality);
						//}
					//} break; 
				//} case 6 : {
					//switch(setting) {
						//case 1: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 0);} break;
						//case 2: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 1);} break;
						//case 3: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 2);} break;
						//case 4: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 3);} break;
						//case 5: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 4);} break;
						//case 6: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 5);} break;
						//case 7: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 6);} break;
						//case 8: { savePMTHitsTiming(&runMCS1, &runMCS2, &time2d, &asym2d, 7);} break;
						//default:{
							//printf("Invalid setting %d for function %d!\n",setting, functionality);
						//}
					//} break;
				//} default :
					//printf("Invalid Functionality Specified!\n");
			//}
						
			//fflush(stdout);
							
		//} else {
			//printf("Warning! Run %05d does not exist! \n", runNo);
		//}
		
		////delete *runBody;
	//}
	////}
	//// Saving externals...
	
	
	//// Temporary doing thing. Going to make this a fucking mess for a bit.
	//switch (functionality) {
		//case 1 : {	} break;
		//case 2 : {
		///*	// Initialize saving names
			//char dirName[256];
			//snprintf(dirName, 256, "%s/%%s", saveLoc);
			//char rootName[256];
			//char tmpName[256];
			//char csvName[256];
						
			//switch(setting) {
				//case 1: {
					//// Create character names
					//snprintf(rootName,256,dirName,"summedDipShort.root");
					//snprintf(csvName,256,dirName,"summedDip.csv");

					//dipSumS.SaveAs(rootName); 
					//FILE* outfile = fopen(csvName,"a");
					//for (int ii = 0; ii < 2099; ii++) {
						//if ((dipSumS.GetBinContent(ii) > 0) && (dipSumS.GetBinContent(ii+1) > 0)) {
							//fprintf(outfile, "%d,%d\n", (int)(dipSumS.GetBinContent(ii)),ii);
						//} else {
							//fprintf(outfile, "%d,%d\n",0,ii);
						//}
					//}
					//printf("Outputting data into .csv form in '%s'\n",csvName);
					//fclose(outfile);				
					//} break;
				//case 2: {
					//char rootName[256];
					//snprintf(rootName,256,dirName,"summedDipLong.root");
					//dipSumL.SaveAs(rootName); 
					//} break;
				//case 10: {
					//// Hard code this stuff
					
					//snprintf(rootName,256,dirName,"fillSumD1.root");
					//fillSumD1.SaveAs(rootName);
					
					//snprintf(rootName,256,dirName,"fillSumD2.root");
					//fillSumD2.SaveAs(rootName);
					
					//snprintf(rootName,256,dirName,"fillSumOl.root");
					//fillSumOl.SaveAs(rootName);
					
					//snprintf(rootName,256,dirName,"fillSumBa.root");
					//fillSumBa.SaveAs(rootName);
					
					//snprintf(rootName,256,dirName,"fillSumSp.root");
					//fillSumSp.SaveAs(rootName);
					//} break;
				//case 4: {
					//// Create character names
					//if (dipSumS.GetBinContent(0) > 0) {
						//printf("Saving Unload\n");
						
						//snprintf(tmpName,256,dirName,"summedRemove.root%d");
						//snprintf(rootName,256,tmpName,jobid);
						//snprintf(tmpName,256,dirName,"oneRemove.csv%d");
						//snprintf(csvName,256,tmpName,jobid);
						//dipSumS.SaveAs(rootName);
						
						//FILE* outfile = fopen(csvName,"w");
						//for (int ii = 0; ii < 3100; ii++) {
							//fprintf(outfile,"%d,%d\n",ii,(int)dipSumS.GetBinContent(ii));
						//}
						//fclose(outfile);
					//}
					//if (dayBkg.GetBinContent(0) > 0) {
						//printf("Saving Background\n");
						//snprintf(tmpName,256,dirName,"summedRemoveBkg.root%d");
						//snprintf(rootName,256,tmpName,jobid);
						//snprintf(tmpName,256,dirName,"oneRemoveBkg.csv%d");
						//snprintf(csvName,256,tmpName,jobid);
						//dayBkg.SaveAs(rootName);
						//FILE* outfile2 = fopen(csvName,"w");
						//for (int ii = 0; ii < 10000; ii++) {
							//fprintf(outfile2,"%d,%d\n",ii,(int)dayBkg.GetBinContent(ii));
						//}
						//fclose(outfile2);
					//}
					//} break;
				//case 5: {
					
					//} break;
				//case 6: {
					//snprintf(rootName,256,dirName,"acSumAll");
					//acSumAll.SaveAs(rootName);
					//} break;
				//case 7: {
					//// Save the root files
					//snprintf(tmpName,256,dirName,"longHold1.root%d");
					//snprintf(rootName,256,tmpName,jobid);
					//longPMT1.SaveAs(rootName);
					//snprintf(tmpName,256,dirName,"longHold2.root%d");
					//snprintf(rootName,256,tmpName,jobid);
					//longPMT2.SaveAs(rootName);
					//snprintf(tmpName,256,dirName,"longHoldC.root%d");
					//snprintf(rootName,256,tmpName,jobid);
					//longCoinc.SaveAs(rootName);
					
					//// Save an output .csv file
					//snprintf(tmpName,256,dirName,"holdingBkg.csv%d");
					//snprintf(csvName,256,tmpName,jobid);
					//FILE* outfile = fopen(csvName, "w");
					//for (int ii = 0; ii < 1550; ii++) {
						//fprintf(outfile,"\n%d",ii);
						//fprintf(outfile,",%f",longPMT1.GetBinContent(ii));
						//fprintf(outfile,",%f",longPMT2.GetBinContent(ii));
						//fprintf(outfile,",%f",longCoinc.GetBinContent(ii));
					//}					
					//fclose(outfile);
					//} break;
				//case 8: {
					//snprintf(rootName,256,dirName,"bkgSum1.root");
					//bkgSum1.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"bkgSum2.root");
					//bkgSum2.SaveAs(rootName);
					//} break;
				//case 9: {
					//snprintf(rootName,256,dirName,"dipSumACS.root");
					//dipSumACS.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"dipSumS.root");
					//dipSumS.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"dipSumACL.root");
					//dipSumACL.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"dipSumL.root");
					//dipSumL.SaveAs(rootName);
					//}break;
				//case 11: {
					//// Save the root files
					//snprintf(tmpName,256,dirName,"totalPMT1.root%d");
					//snprintf(rootName,256,tmpName,jobid);
					//totalPMT1.SaveAs(rootName);
					//snprintf(tmpName,256,dirName,"totalPMT2.root%d");
					//snprintf(rootName,256,tmpName,jobid);
					//totalPMT2.SaveAs(rootName);
					//snprintf(tmpName,256,dirName,"totalPMTC.root%d");
					//snprintf(rootName,256,tmpName,jobid);
					//totalCoin.SaveAs(rootName);
					
					//// Save an output .csv file
					//snprintf(tmpName,256,dirName,"totalByPMT.csv%d");
					//snprintf(csvName,256,tmpName,jobid);
					//FILE* outfile = fopen(csvName, "w");
					//for (int ii = 0; ii < 1000; ii++) {
						//fprintf(outfile,"\n%d",ii);
						//fprintf(outfile,",%f",totalPMT1.GetBinContent(ii));
						//fprintf(outfile,",%f",totalPMT2.GetBinContent(ii));
						//fprintf(outfile,",%f",totalCoin.GetBinContent(ii));
					//}					
					//fclose(outfile);	
				//}break;
				//default: {printf("Invalid setting %d for function %d!\n",setting, functionality); } break;
			 
			//}*/
		//} break;
		//case 3 : {	} break;
		//case 4 : {
			//// Initialize saving names
			//char dirName[256];
			//snprintf(dirName, 256, "%s/%%s", saveLoc);
			//char rootName[256];
			//char csvName[256];
			//char tmpName[256];	
			//switch(setting) {
				//case 1: { 
					//// Save histograms (ROOT)
					////std::string tmpName;
					//snprintf(rootName,256,dirName,"pTiming1Bkg.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming1.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming2Bkg.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming2.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming3Bkg.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming3.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming4Bkg.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming4.SaveAs(tmpName);
					//phsHits.Draw("COLZ");
					//snprintf(rootName,256,dirName,"phsHitsBkg.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//phsHits.SaveAs(tmpName);
					
					//int binNo = 0;
					
					//// Convert histograms into .csv file
					
					//snprintf(csvName,256,dirName,"phsHitsBkg.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfile = fopen(tmpName, "w");
					//for (int ii = 1; ii <= phsTB; ii++) {
						//for (int jj = 1; jj <= phsTB; jj++) {
							//binNo = phsHits.GetBin(ii,jj);
							//fprintf(outfile,"%d,%d:%d\n",(int)phsHits.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfile);
					//// Save timing histos as csv
					//snprintf(csvName,256,dirName,"pTiming1Bkg.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile1 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming2Bkg.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile2 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming3Bkg.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile3 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming4Bkg.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile4 = fopen(tmpName,"w");
					//for (int ii = 1; ii <= nTB; ii++) {
						//fprintf(pOutfile1,"%d,%d\n",(int)pTiming1.GetBinContent(ii),ii);
						//fprintf(pOutfile2,"%d,%d\n",(int)pTiming2.GetBinContent(ii),ii);
						//fprintf(pOutfile3,"%d,%d\n",(int)pTiming3.GetBinContent(ii),ii);
						//fprintf(pOutfile4,"%d,%d\n",(int)pTiming4.GetBinContent(ii),ii);
					//}
					//fclose(pOutfile1);
					//fclose(pOutfile2);
					//fclose(pOutfile3);
					//fclose(pOutfile4);
					
					//} break;
				//case 2: {
					//// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"pTiming1.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming1.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming2.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming2.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming3.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming3.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming4.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming4.SaveAs(tmpName);
					//phsHits.Draw("COLZ");
					//snprintf(rootName,256,dirName,"phsHits.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//phsHits.SaveAs(tmpName);
					
					//int binNo = 0;
					
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"phsHits.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfile = fopen(tmpName, "w");
					//for (int ii = 1; ii <= phsTB; ii++) {
						//for (int jj = 1; jj <= phsTB; jj++) {
							//binNo = phsHits.GetBin(ii,jj);
							//fprintf(outfile,"%d,%d:%d\n",(int)phsHits.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfile);
					//// Save timing histos as csv
					//snprintf(csvName,256,dirName,"pTiming1.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile1 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming2.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile2 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming3.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile3 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming4.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile4 = fopen(tmpName,"w");
					//for (int ii = 1; ii <= nTB; ii++) {
						//fprintf(pOutfile1,"%d,%d\n",(int)pTiming1.GetBinContent(ii),ii);
						//fprintf(pOutfile2,"%d,%d\n",(int)pTiming2.GetBinContent(ii),ii);
						//fprintf(pOutfile3,"%d,%d\n",(int)pTiming3.GetBinContent(ii),ii);
						//fprintf(pOutfile4,"%d,%d\n",(int)pTiming4.GetBinContent(ii),ii);
					//}
					//fclose(pOutfile1);
					//fclose(pOutfile2);
					//fclose(pOutfile3);
					//fclose(pOutfile4);
					//} break;
				//case 3: {
					//// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"pTiming1HoldS.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming1.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming2HoldS.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming2.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming3HoldS.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming3.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming4HoldS.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming4.SaveAs(tmpName);
					//phsHits.Draw("COLZ");
					//snprintf(rootName,256,dirName,"phsHitsHoldS.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//phsHits.SaveAs(tmpName);
					
					//int binNo = 0;
					
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"phsHitsHoldS.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfile = fopen(tmpName, "w");
					//for (int ii = 1; ii <= phsTB; ii++) {
						//for (int jj = 1; jj <= phsTB; jj++) {
							//binNo = phsHits.GetBin(ii,jj);
							//fprintf(outfile,"%d,%d:%d\n",(int)phsHits.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfile);
					//// Save timing histos as csv
					//snprintf(csvName,256,dirName,"pTiming1HoldS.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile1 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming2HoldS.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile2 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming3HoldS.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile3 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming4HoldS.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile4 = fopen(tmpName,"w");
					//for (int ii = 1; ii <= nTB; ii++) {
						//fprintf(pOutfile1,"%d,%d\n",(int)pTiming1.GetBinContent(ii),ii);
						//fprintf(pOutfile2,"%d,%d\n",(int)pTiming2.GetBinContent(ii),ii);
						//fprintf(pOutfile3,"%d,%d\n",(int)pTiming3.GetBinContent(ii),ii);
						//fprintf(pOutfile4,"%d,%d\n",(int)pTiming4.GetBinContent(ii),ii);
					//}
					//fclose(pOutfile1);
					//fclose(pOutfile2);
					//fclose(pOutfile3);
					//fclose(pOutfile4);
					//} break;
				//case 4 : { 
				//// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"pTiming1Unl.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming1.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming2Unl.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming2.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming3Unl.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming3.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"pTiming4Unl.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//pTiming4.SaveAs(tmpName);
					//phsHits.Draw("COLZ");
					//snprintf(rootName,256,dirName,"phsHitsUnl.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//phsHits.SaveAs(tmpName);
					
					//int binNo = 0;
					
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"phsHitsUnl.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfile = fopen(tmpName, "w");
					//for (int ii = 1; ii <= phsTB; ii++) {
						//for (int jj = 1; jj <= phsTB; jj++) {
							//binNo = phsHits.GetBin(ii,jj);
							//fprintf(outfile,"%d,%d:%d\n",(int)phsHits.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfile);
					//// Save timing histos as csv
					//snprintf(csvName,256,dirName,"pTiming1Unl.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile1 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming2Unl.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile2 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming3Unl.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile3 = fopen(tmpName,"w");
					//snprintf(csvName,256,dirName,"pTiming4Unl.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* pOutfile4 = fopen(tmpName,"w");
					//for (int ii = 1; ii <= nTB; ii++) {
						//fprintf(pOutfile1,"%d,%d\n",(int)pTiming1.GetBinContent(ii),ii);
						//fprintf(pOutfile2,"%d,%d\n",(int)pTiming2.GetBinContent(ii),ii);
						//fprintf(pOutfile3,"%d,%d\n",(int)pTiming3.GetBinContent(ii),ii);
						//fprintf(pOutfile4,"%d,%d\n",(int)pTiming4.GetBinContent(ii),ii);
					//}
					//fclose(pOutfile1);
					//fclose(pOutfile2);
					//fclose(pOutfile3);
					//fclose(pOutfile4);
					//} break;
				//case 5: {
					//snprintf(rootName,256,dirName,"pmt1DT.root"+jobid);
					//pmt1DT.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"pmt2DT.root"+jobid);
					//pmt2DT.SaveAs(rootName);		
					//} break;
				//case 6: {
					//snprintf(rootName,256,dirName,"pmt1DT.root"+jobid);
					//pmt1DT.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"pmt2DT.root"+jobid);
					//pmt2DT.SaveAs(rootName);					
				//} break;
				//case 7: {// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"pT2D1_bkg.root"+jobid);
					//pT2D1.Draw("COLZ");
					//pT2D1.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"pT2D2_bkg.root"+jobid);
					//pT2D2.Draw("COLZ");
					//pT2D2.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"pT2D3_bkg.root"+jobid);
					//pT2D3.Draw("COLZ");
					//pT2D3.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"pT2D4_bkg.root"+jobid);
					//pT2D4.Draw("COLZ");
					//pT2D4.SaveAs(rootName);
					
					//int binNo = 0;
					
					//// Save timing histos as csv
					//snprintf(csvName,256,dirName,"pT2D1_bkg.csv"+jobid);
					//FILE* pOut2DB1 = fopen(csvName,"w");
					//snprintf(csvName,256,dirName,"pT2D2_bkg.csv"+jobid);
					//FILE* pOut2DB2 = fopen(csvName,"w");
					//snprintf(csvName,256,dirName,"pT2D3_bkg.csv"+jobid);
					//FILE* pOut2DB3 = fopen(csvName,"w");
					//snprintf(csvName,256,dirName,"pT2D4_bkg.csv"+jobid);
					//FILE* pOut2DB4 = fopen(csvName,"w");
					//for (int ii = 1; ii <= nTB; ii++) {
						//for (int jj = 1; jj <= phsTB; jj++) {
							//binNo = pT2D1.GetBin(ii,jj);
							//fprintf(pOut2DB1,"%d,%d,%d\n",(int)pT2D1.GetBinContent(binNo),ii,jj);
							//fprintf(pOut2DB2,"%d,%d,%d\n",(int)pT2D2.GetBinContent(binNo),ii,jj);
							//fprintf(pOut2DB3,"%d,%d,%d\n",(int)pT2D3.GetBinContent(binNo),ii,jj);
							//fprintf(pOut2DB4,"%d,%d,%d\n",(int)pT2D4.GetBinContent(binNo),ii,jj);
						//}
					//}
					//fclose(pOut2DB1);
					//fclose(pOut2DB2);
					//fclose(pOut2DB3);
					//fclose(pOut2DB4);
					//} break;
				//case 8: {
					//snprintf(rootName,256,dirName,"pT2D1.root"+jobid);
					//pT2D1.Draw("COLZ");
					//pT2D1.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"pT2D2.root"+jobid);
					//pT2D2.Draw("COLZ");
					//pT2D2.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"pT2D3.root"+jobid);
					//pT2D3.Draw("COLZ");
					//pT2D3.SaveAs(rootName);
					//snprintf(rootName,256,dirName,"pT2D4.root"+jobid);
					//pT2D4.Draw("COLZ");
					//pT2D4.SaveAs(rootName);
					
					//int binNo = 0;
					
					//// Save timing histos as csv
					//snprintf(csvName,256,dirName,"pT2D1.csv"+jobid);
					//FILE* pOut2D1 = fopen(csvName,"w");
					//snprintf(csvName,256,dirName,"pT2D2.csv"+jobid);
					//FILE* pOut2D2 = fopen(csvName,"w");
					//snprintf(csvName,256,dirName,"pT2D3.csv"+jobid);
					//FILE* pOut2D3 = fopen(csvName,"w");
					//snprintf(csvName,256,dirName,"pT2D4.csv"+jobid);
					//FILE* pOut2D4 = fopen(csvName,"w");
					//for (int ii = 1; ii <= nTB; ii++) {
						//for (int jj = 1; jj <= phsTB; jj++) {
							//binNo = pT2D1.GetBin(ii,jj);
							//fprintf(pOut2D1,"%d,%d,%d\n",(int)pT2D1.GetBinContent(binNo),ii,jj);
							//fprintf(pOut2D2,"%d,%d,%d\n",(int)pT2D2.GetBinContent(binNo),ii,jj);
							//fprintf(pOut2D3,"%d,%d,%d\n",(int)pT2D3.GetBinContent(binNo),ii,jj);
							//fprintf(pOut2D4,"%d,%d,%d\n",(int)pT2D4.GetBinContent(binNo),ii,jj);
						//}
					//}
					//fclose(pOut2D1);
					//fclose(pOut2D2);
					//fclose(pOut2D3);
					//fclose(pOut2D4);
					//} break;
				//case 9: {
					//snprintf(csvName,256,"fillDT_%d.root",coincMode);
					//snprintf(rootName,256, dirName,csvName);
					//pmt1DT.SaveAs(rootName);
				//} break;
				//default : {printf("Invalid setting %d for function %d!\n",setting, functionality); } break;
			//}
		//} break;	
		//case 5 : {
			//// Initialize saving names
			//char dirName[256];
			//snprintf(dirName, 256, "%s/pmtHits/%%s", saveLoc);
			//char rootName[256];
			//switch(setting) {
				//case 1: { 
						//// Save histograms (ROOT)
						//snprintf(rootName, 256, dirName,"lTMax.root");
						//lTMax.SaveAs(rootName);
						//snprintf(rootName, 256, dirName,"lNPE.root");
						//lNPE.SaveAs(rootName);
				//} break;
				//case 2: {
					
					//snprintf(rootName, 256, dirName,"numPh.root");
					//numPh.SaveAs(rootName);
					
					//snprintf(rootName, 256, dirName,"coinLen.root");
					//coinLen.SaveAs(rootName);
										
					//snprintf(rootName, 256, dirName,"phByLen.root");
					//phByLen.SaveAs(rootName);
				//} break;
			//default : {printf("Invalid setting %d for function %d!\n",setting, functionality); } break;
			//}
		//} break;
		//case 6 : { 
			//// Initialize saving names
			//char dirName[256];
			//snprintf(dirName, 256, "%s/%%s", saveLoc);
			//char rootName[256];
			//char csvName[256];
			//char tmpName[256];	
			//switch(setting) {
				//case 1: { 
					//// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"hitsByLenBkg38.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//time2d.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"asymByLenBkg38.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//asym2d.SaveAs(tmpName);
					
					//int binNo = 0;
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"hitsByLenBkg38.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileH = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = time2d.GetBin(ii,jj);
							//fprintf(outfileH,"%d,%d:%d\n",(int)time2d.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfileH);
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"asymByLenBkg38.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileA = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = asym2d.GetBin(ii,jj);
							//fprintf(outfileA,"%d,%d:%d\n",(int)asym2d.GetBinContent(binNo),ii,jj);
						//}
					//}	
					//fclose(outfileA);
				//} break;
				//case 2: {
					//// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"hitsByLenBkg49.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//time2d.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"asymByLenBkg49.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//asym2d.SaveAs(tmpName);
					
					//int binNo = 0;
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"hitsByLenBkg49.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileH = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = time2d.GetBin(ii,jj);
							//fprintf(outfileH,"%d,%d:%d\n",(int)time2d.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfileH);
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"asymByLenBkg49.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileA = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = asym2d.GetBin(ii,jj);
							//fprintf(outfileA,"%d,%d:%d\n",(int)asym2d.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfileA);
				//} break;
				//case 3: {
					//// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"hitsByLenBkg25.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//time2d.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"asymByLenBkg25.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//asym2d.SaveAs(tmpName);
					
					//int binNo = 0;
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"hitsByLenBkg25.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileH = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = time2d.GetBin(ii,jj);
							//fprintf(outfileH,"%d,%d:%d\n",(int)time2d.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfileH);
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"asymByLenBkg25.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileA = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = asym2d.GetBin(ii,jj);
							//fprintf(outfileA,"%d,%d:%d\n",(int)asym2d.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfileA);
				//} break;
				//case 4 : { 
					//// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"hitsByLenBkg1.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//time2d.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"asymByLenBkg1.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//asym2d.SaveAs(tmpName);
					
					//int binNo = 0;
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"hitsByLenBkg1.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileH = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = time2d.GetBin(ii,jj);
							//fprintf(outfileH,"%d,%d:%d\n",(int)time2d.GetBinContent(binNo),ii,jj);
						//}
					//}	
					//fclose(outfileH);
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"asymByLenBkg1.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileA = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = asym2d.GetBin(ii,jj);
							//fprintf(outfileA,"%d,%d:%d\n",(int)asym2d.GetBinContent(binNo),ii,jj);
						//}
					//}
					//fclose(outfileA);
				//} break;
				//case 5: { 
					//// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"hitsByLenBkg.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//time2d.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"asymByLenBkg.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//asym2d.SaveAs(tmpName);
					
					//int binNo = 0;
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"hitsByLenBkg.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileH = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = time2d.GetBin(ii,jj);
							//fprintf(outfileH,"%d,%d:%d\n",(int)time2d.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfileH);
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"asymByLenBkg.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileA = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = asym2d.GetBin(ii,jj);
							//fprintf(outfileA,"%d,%d:%d\n",(int)asym2d.GetBinContent(binNo),ii,jj);
						//}
					//}	
					//fclose(outfileA);
				//} break;
				//case 6: {
					//// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"hitsByLen.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//time2d.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"asymByLen.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//asym2d.SaveAs(tmpName);
					
					//int binNo = 0;
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"hitsByLen.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileH = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = time2d.GetBin(ii,jj);
							//fprintf(outfileH,"%d,%d:%d\n",(int)time2d.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfileH);
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"asymByLen.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileA = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = asym2d.GetBin(ii,jj);
							//fprintf(outfileA,"%d,%d:%d\n",(int)asym2d.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfileA);
				//} break;
				//case 7 : {
					//// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"hitsByLenHold.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//time2d.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"asymByLenHold.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//asym2d.SaveAs(tmpName);
					
					//int binNo = 0;
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"hitsByLenHold.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileH = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = time2d.GetBin(ii,jj);
							//fprintf(outfileH,"%d,%d:%d\n",(int)time2d.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfileH);
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"asymByLenHold.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileA = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = asym2d.GetBin(ii,jj);
							//fprintf(outfileA,"%d,%d:%d\n",(int)asym2d.GetBinContent(binNo),ii,jj);
						//}
					//}					
					//fclose(outfileA);
				//} break;
				//case 8 : { 
					//// Save histograms (ROOT)
					//snprintf(rootName,256,dirName,"hitsByLenUnl.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//time2d.SaveAs(tmpName);
					//snprintf(rootName,256,dirName,"asymByLenUnl.root%d");
					//snprintf(tmpName,256,rootName,jobid);
					//asym2d.SaveAs(tmpName);
					
					//int binNo = 0;
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"hitsByLenUnl.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileH = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = time2d.GetBin(ii,jj);
							//fprintf(outfileH,"%d,%d:%d\n",(int)time2d.GetBinContent(binNo),ii,jj);
						//}
					//}	
					//fclose(outfileH);
					//// Convert histograms into .csv file
					//snprintf(csvName,256,dirName,"asymByLenUnl.csv%d");
					//snprintf(tmpName,256,csvName,jobid);
					//FILE* outfileA = fopen(tmpName, "w");
					//for (int ii = 1; ii <= bin2d; ii++) {
						//for (int jj = 1; jj <= bin2d; jj++) {
							//binNo = asym2d.GetBin(ii,jj);
							//fprintf(outfileA,"%d,%d:%d\n",(int)asym2d.GetBinContent(binNo),ii,jj);
						//}
					//}
					//fclose(outfileA);
				//} break;
			//} break;
		//} break; 			
		//default : 
			//printf("Not saving external ROOT files!\n");
	//}	
	//return 0;
}

		/*// Hard coding in the 2017/2018 year split at 9519
		// Also hard coding in run 3412 as the beginning of 2017.
		// 9519-9538 are in both directories BUT only in 2018 runlogdump.
		if ((3412 <= runNo) && (runNo < 9519)) {
			snprintf(fName, 256, "%s/processed_output_%%05d.root", fileLoc17);
		} else if ((9519 <= runNo) && (runNo < 14732)) {
			snprintf(fName, 256, "%s/processed_output_%%05d.root", fileLoc18);
		} else if((14732 <= runNo) && (runNo < 17999)) {
			snprintf(fName, 256, "%s/processed_output_%%05d.root", fileLoc19);
		} else {
			printf("Warning! Run %05d is not from calendar year 2017 or 2018 !\n", runNo);
			continue;
		}
				
		std::string runBody = fName;
		sprintf(fName, runBody.c_str(), runNo);
		std::ifstream infile(fName);
		
		if (infile.good()) {
			Run runMCS1(coincWindow, peSumWindow, peSum, runNo, coincMode, fName, "tmcs_0");
			Run runMCS2(coincWindow, peSumWindow, peSum, runNo, coincMode, fName, "tmcs_1");
			EMS runEMS(runNo, fName, "tems");
			
			// MCS2 changes by 10 at run 9020:
			if (runNo >= 9020) {
				runMCS2.setMCSOff(10);
			}
	*/
/*			switch(setting) {
				case 1: {normNByDip(&runMCS1, &runMCS2, &runEMS, traceLoc);     } break;
				case 2: {normNByDipSing(&runMCS1, &runMCS2, &runEMS, traceLoc); } break;
				case 3: {normNByDipCoinc(&runMCS1, &runMCS2, &runEMS, traceLoc);} break;
				case 4: {extractObservables(&runMCS1, &runMCS2, traceLoc);} break;
				case 5: {activeCleaner(&runMCS1,&runMCS2,&runEMS,traceLoc);} break;
				case 6: {fitBkgUnload(&runMCS1,&runMCS2,&runEMS);} break;
				case 7: {normNByDipDet(&runMCS1,&runMCS2,&runEMS,traceLoc);} break;
				default:{
					printf("Invalid setting %d for function %d!\n",setting, functionality);
				} 
			} break;
		}
	
}*/











//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* Removed (commented) code for cleanliness: (Should eventually decide whether 
 * to ultimately remove this or put it back in.) (Line numbers might be inaccurate.)
 * 		TH1D ttne("ttne", "ttne", 1000,0,1000);
		TH1D ttpe("ttpe", "ttpe", 1000,0,1000);
			
	//if(strstr(query.c_str(), "SELECT")) {
	//	DBHandler hand(query.c_str(), coincWindow, peSumWindow, peSum, coincMode);
	//} else {
	 * 47(between double + std::vector): printf("fillEnd: %f\n", fillEnd);
	 * 53(after std::for_each): //spHits.insert(spHits.end(), spCts.begin(), spCts.end());
	 * 63(before TH1D singlebkg)
		auto dipSummer = [&sCts, &lCts](Run* run) {
		double firstDip = run->getTagBitEvt(1<<9, 175, 0);
		std::vector<input_t> dCts = run->getCoincCounts(
			[firstDip](input_t x)->input_t{x.realtime -= firstDip; return x;},
			[firstDip](input_t x)->bool{return true;}
			//[firstDip](input_t x)->bool{return x.realtime > firstDip;}
		);
		std::vector<input_t> dCts = run->getCoincCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return true;}
		);
		//if(firstDip < 1000) {
		 if(dCts.empty()) {
			return;
		 }
		 if(dCts.back().realtime > 0 && dCts.back().realtime < 100000) {
			printf("Added Short\n");
			std::for_each(dCts.begin(), dCts.end(), [&sCts](input_t x){sCts.push_back(x.realtime);});
		 } 
		 else {
			printf("Added Long\n");
			std::for_each(dCts.begin(), dCts.end(), [&lCts](input_t x){lCts.push_back(x.realtime);});
			}
		};
	 * 80(after std::for_each(dCts.begin(), dCts.end(),...):
	    auto ttneAC = [&ttne, &ttpe](Run* mcs1, Run* mcs2) {
		struct coinc {
			double t;
			int ch;
		};
		//double firstDip = mcs1->getTagBitEvt(1<<9, 175, 0);
		double firstDip = 0.0;
		std::vector<input_t> dCts = mcs1->getCoincCounts(
			[firstDip](input_t x)->input_t{x.realtime -= firstDip; return x;},
			[firstDip](input_t x)->bool{return x.realtime > firstDip && x.realtime < firstDip + 50;}
		);
		std::vector<input_t> acCts = mcs2->getCoincCounts(
			[firstDip](input_t x)->input_t{x.realtime -= firstDip; return x;},
			[firstDip](input_t x)->bool{return x.realtime > firstDip && x.realtime < firstDip + 50;}
		);
		std::vector<coinc> allCts;
		if(dCts.empty() || acCts.empty()) {
			printf("Empty! dCts: %ld aCts: %ld\n", dCts.size(), acCts.size());
			return;
		}
		std::for_each(dCts.begin(), dCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 1}; allCts.push_back(y);});
		std::for_each(acCts.begin(), acCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 2}; allCts.push_back(y);});
		std::sort(allCts.begin(), allCts.end(), [](coinc x, coinc y)->bool{return (x.t < y.t);});
		for(auto it = allCts.begin() ; it < allCts.end(); it++) {
			if(it->ch != 2) { continue; }
			for(auto nIt = it; nIt < allCts.end(); nIt++) {
				if(nIt->ch == 1) {
					ttne.Fill(((nIt->t - it->t) / NANOSECOND) / 1000.0);
					break;
				}
			}
			for(auto prevIt = it; prevIt > allCts.begin(); prevIt--) {
				if(prevIt->ch == 1) {
					ttpe.Fill(((it->t - prevIt->t) / NANOSECOND) / 1000.0);
					break;
				}
			}
		}
	}
	* /*auto coincACSummer = [&sCts, &lCts](Run* mcs1, Run* mcs2) {
		struct coinc {
			double t;
			int ch;
		};
		double firstDip = mcs1->getTagBitEvt(1<<9, 175, 0);
		//double firstDip = 0.0;
		std::vector<input_t> dCts = mcs1->getCoincCounts(
			[firstDip](input_t x)->input_t{x.realtime -= firstDip; return x;},
			[firstDip](input_t x)->bool{return x.realtime > firstDip && x.realtime < firstDip + 50;}
		);
		std::vector<input_t> acCts = mcs2->getCoincCounts(
			[firstDip](input_t x)->input_t{x.realtime -= firstDip; return x;},
			[firstDip](input_t x)->bool{return x.realtime > firstDip && x.realtime < firstDip + 50;}
		);
		std::vector<input_t> antiCoincAC;
		std::vector<coinc> allCts;
		if(dCts.empty() || acCts.empty()) {
			printf("Empty! dCts: %ld aCts: %ld\n", dCts.size(), acCts.size());
			return;
		}
		std::for_each(dCts.begin(), dCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 1}; allCts.push_back(y);});
		std::for_each(acCts.begin(), acCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 2}; allCts.push_back(y);});
		std::sort(allCts.begin(), allCts.end(), [](coinc x, coinc y)->bool{return (x.t < y.t);});
		for(auto it = allCts.begin() ; it < allCts.end(); it++) {
			if(it->ch != 2) { continue; }
			for(auto nIt = it + 1; nIt < allCts.end(); nIt++) {
				if(nIt->ch == 1 && (nIt->t - it->t) < 60000 * NANOSECOND) {
					break;
				}
				else if(nIt->ch == 1) {
					input_t event = {0, it->t, 2, 0};
					antiCoincAC.push_back(event);
					break;
				}
			}
		}
		if(antiCoincAC.empty()) {
			return;
		}
		if(antiCoincAC.back().realtime > 500 && antiCoincAC.back().realtime < 1000) {
			printf("Added Short\n");
			std::for_each(antiCoincAC.begin(), antiCoincAC.end(), [&sCts](input_t x){sCts.push_back(x.realtime);});
			std::for_each(dCts.begin(), dCts.end(), [&lCts](input_t x){lCts.push_back(x.realtime);});
		}
		else {
			printf("Added Long\n");
			std::for_each(antiCoincAC.begin(), antiCoincAC.end(), [&lCts](input_t x){lCts.push_back(x.realtime);});
			std::for_each(dCts.begin(), dCts.end(), [&sCts](input_t x){sCts.push_back(x.realtime);});
		}
	};*/
	/*auto coincACIntegral = [](Run* mcs1, Run* mcs2) {
		int numTail = 0;
		int numHold = 0;
		int numCount = 0;
		struct coinc {
			double t;
			int ch;
		};
		double fillEnd = mcs1->getTagBitEvt(1<<9, 50, 0);
		double firstDip = mcs1->getTagBitEvt(1<<9, 175, 0);
		double secondDip = mcs1->getTagBitEvt(1<<9, firstDip + 10, 0);
		std::vector<input_t> dCts = mcs1->getCoincCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return true;}
		);
		std::vector<input_t> acCts = mcs2->getCoincCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return true;}
		);
		std::vector<coinc> allCts;
		if(dCts.empty() || acCts.empty()) {
			printf("Empty! dCts: %ld aCts: %ld\n", dCts.size(), acCts.size());
			return;
		}
		std::for_each(dCts.begin(), dCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 1}; allCts.push_back(y);});
		std::for_each(acCts.begin(), acCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 2}; allCts.push_back(y);});
		std::sort(allCts.begin(), allCts.end(), [](coinc x, coinc y)->bool{return (x.t < y.t);});
		for(auto it = allCts.begin() ; it < allCts.end(); it++) {
			if(it->ch != 2) { continue; }
			for(auto nIt = it + 1; nIt < allCts.end(); nIt++) {
				if(nIt->ch == 1 && (nIt->t - it->t) < 60000 * NANOSECOND) {
					break;
				}
				else if(nIt->ch == 1) {
					if(nIt->t > fillEnd + 100.0 && nIt->t < fillEnd + 200.0) { numTail += 1; }
					if(nIt->t > fillEnd + 310.0 && nIt->t < firstDip) { numHold += 1; }
					if(nIt->t > firstDip && nIt->t < secondDip) { numCount += 1; }
					break;
				}
			}
		}
		printf("Data - %d,%f,%d,%f,%d,%f\n", numTail, 100.0, numHold, firstDip - (fillEnd + 310.0), numCount, secondDip - firstDip);
	};*/
	/*TH1D scoreA("scoreA", "scoreA", 10000, 0, 10000);
	TH1D scoreB("scoreB", "scoreB", 10000, 0, 10000);
	auto noiseFinder = [&scoreA, &scoreB](Run* run) {
		std::vector<input_t> ctsA = run->getCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return x.ch == 1;}
		);
		std::vector<input_t> ctsB = run->getCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return x.ch == 2;}
		);
		for(auto it = ctsA.begin(); it < ctsA.end()-2; it++) {
			double scoreSum = 0.0;
			for(int nPeriod = 1; nPeriod <= 10; nPeriod++) {
				double prevDelta = ((it->realtime + nPeriod*(1.0/20000.0)) - (it+1)->realtime)/(1.0E-6);
				for(auto cIt = it+2; cIt < ctsA.end(); cIt++) {
					double delta = ((it->realtime + nPeriod*(1.0/20000.0)) - cIt->realtime)/(1.0E-6);
					if(fabs(delta) > fabs(prevDelta)) {
						break;
					}
					prevDelta = delta;
				}
				if(fabs(prevDelta) > 0.0) { scoreSum += 1.0/(prevDelta*prevDelta); }
				if(1.0/(prevDelta*prevDelta) > 6100 && 1.0/(prevDelta*prevDelta) < 6120) {
					printf("%.18f\n", prevDelta);
				}
			}
			scoreA.Fill(scoreSum);
		}
		for(auto it = ctsB.begin(); it < ctsB.end()-2; it++) {
			double scoreSum = 0.0;
			for(int nPeriod = 1; nPeriod <= 10; nPeriod++) {
				double prevDelta = ((it->realtime + nPeriod*(1.0/20000.0)) - (it+1)->realtime)/(1.0E-6);
				for(auto cIt = it+2; cIt < ctsB.end(); cIt++) {
					double delta = ((it->realtime + nPeriod*(1.0/20000.0)) - cIt->realtime)/(1.0E-6);
					if(fabs(delta) > fabs(prevDelta)) {
						break;
					}
					prevDelta = delta;
				}
				if(fabs(prevDelta) > 0.0) { scoreSum += 1.0/(prevDelta*prevDelta); }
			}
			scoreB.Fill(scoreSum);
		}
	};
	; */
	// These objects exist inside our main()
	// They contain different counting vectors/histograms
	// They're here (and not in Functions.cpp) because these are total sums.
	
	//------------------------------------------------------------------
	// Filling sum histograms:
	// Dagger PMT1
	/*auto fillSummerD1 = [&fillSumD1](Run* mcs1, double fillEnd) {
		
		// Make sure we	have a real fillEnd
		if(fillEnd <= 0.0) { 
			return; 
		}
		printf("FillEnd: %f: \n", fillEnd);
		
		// Load the counts from mcs1
		std::vector<input_t> d1Cts = mcs1->getCounts(
			[fillEnd](input_t x)->input_t{input_t y = x; /*y.realtime -= fillEnd; return y;},
			[fillEnd](input_t x)->bool{return x.ch == 1 && x.realtime < fillEnd;}
		);
		if(d1Cts.size() < 100 || d1Cts.empty()) {
			return;
		}

		// Add to our histogram
		printf("Added - %lu to dagger PMT1 \n", d1Cts.size());
		std::for_each(d1Cts.begin(), d1Cts.end(), [&fillSumD1](input_t x){fillSumD1.Fill(x.realtime);});
		
	};
	
	// Dagger PMT2
	auto fillSummerD2 = [&fillSumD2](Run* mcs1, double fillEnd) {
		
		// Make sure we	have a real fillEnd
		if(fillEnd <= 0.0) { 
			return; 
		}
		
		// Load the counts from mcs1
		std::vector<input_t> d2Cts = mcs1->getCounts(
			[fillEnd](input_t x)->input_t{input_t y = x; /*y.realtime -= fillEnd; return y;},
			[fillEnd](input_t x)->bool{return x.ch == 2 && x.realtime < fillEnd;}
		);
		if(d2Cts.size() < 100 || d2Cts.empty()) {
			return;
		}

		// Add to our histogram
		printf("Added - %lu to dagger PMT2 \n", d2Cts.size());
		std::for_each(d2Cts.begin(), d2Cts.end(), [&fillSumD2](input_t x){fillSumD2.Fill(x.realtime);});
		
	};
	
	// Old GV Monitor
	auto fillSummerOl = [&fillSumOl](Run* mcs1, double fillEnd) {
		
		// Make sure we	have a real fillEnd
		if(fillEnd <= 0.0) { 
			return; 
		}
		
		// Load the counts from mcs1
		std::vector<input_t> olCts = mcs1->getCounts(
			[fillEnd](input_t x)->input_t{input_t y = x; /*y.realtime -= fillEnd; return y;},
			[fillEnd](input_t x)->bool{return x.ch == 3 && x.realtime < fillEnd;}
		);
		if(olCts.size() < 100 || olCts.empty()) {
			return;
		}

		// Add to our histogram	
		printf("Added - %lu to old monitor \n", olCts.size());
		std::for_each(olCts.begin(), olCts.end(), [&fillSumOl](input_t x){fillSumOl.Fill(x.realtime);});
	
	};
	
	// Bare (RH?) Monitor
	auto fillSummerBa = [&fillSumBa](Run* mcs1, double fillEnd) {
		
		// Make sure we	have a real fillEnd
		if(fillEnd <= 0.0) { 
			return; 
		}
		
		// Load the counts from mcs1
		std::vector<input_t> baCts = mcs1->getCounts(
			[fillEnd](input_t x)->input_t{input_t y = x; /*y.realtime -= fillEnd; return y;},
			[fillEnd](input_t x)->bool{return x.ch == 4 && x.realtime < fillEnd;}
		);
		if(baCts.size() < 100 || baCts.empty()) {
			return;
		}

		// Add to our histogram
		printf("Added - %lu to bare monitor \n", baCts.size());
		std::for_each(baCts.begin(), baCts.end(), [&fillSumBa](input_t x){fillSumBa.Fill(x.realtime);});
		
	};
	
	// Standpipe Monitor
	auto fillSummerSp = [&fillSumSp](Run* mcs1, double fillEnd) {
		
		// Make sure we	have a real fillEnd
		if(fillEnd <= 0.0) { 
			return; 
		}
		
		// Load the counts from mcs1
		std::vector<input_t> spCts = mcs1->getCounts(
			[fillEnd](input_t x)->input_t{input_t y = x; /*y.realtime -= fillEnd; return y;},
			[fillEnd](input_t x)->bool{return x.ch == 5 && x.realtime < fillEnd;}
		);
		if(spCts.size() < 100 || spCts.empty()) {
			return;
		}

		// Add to our histogram
		printf("Added - %lu to standpipe \n", spCts.size());
		std::for_each(spCts.begin(), spCts.end(), [&fillSumSp](input_t x){fillSumSp.Fill(x.realtime);});
		
	};
	
	//------------------------------------------------------------------
	// Other Summing Histograms
	// Create a vector to count the hits (in standpipe)
	std::vector<double> spHits;
	auto fillSummer = [&spHits](Run* mcs1, Run* mcs2) {
		//double fillEnd = mcs1->getTagBitEvt(1<<0, 290, 0);
		double fillEnd = getFillEnd(mcs1,mcs2);
		std::vector<input_t> spCts = mcs1->getCounts(
			[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
			[fillEnd](input_t x)->bool{return x.ch == 5 && x.realtime < fillEnd;}
		);
		std::for_each(spCts.begin(), spCts.end(), [&spHits](input_t x){spHits.push_back(x.realtime);});
	};
		
	// Fill the background histograms the same way as our vector that created hits
	auto bkgSummer = [&singleBkg](Run* mcs1, Run* mcs2) {
		double firstDip = mcs2->getTagBitEvt(1<<0, 325, 0); // Why is this bit 9?
		std::vector<input_t> dCts = mcs1->getCounts(
			[firstDip](input_t x)->input_t{x.realtime -= (firstDip); return x;},
			[firstDip](input_t x)->bool{return (x.ch==1 || x.ch==2) && x.realtime > (firstDip-500);}
		);
		if(dCts.empty()) {
			return;
		}
		std::for_each(dCts.begin(), dCts.end(), [&singleBkg](input_t x){singleBkg.Fill(x.realtime);});
	};
	
	// Now fill our dip sum histograms
	auto dipSummer = [&dipSumS](Run* mcs1, Run* mcs2) {
		double firstDip = mcs2->getTagBitEvt(1<<0, 325, 1);
		std::vector<double> dagSteps = dagDips(mcs1, mcs2);
		if(dagSteps.size() < 2) { 
			return; 
		}

		// might need to modify position of start/stop for our dips
        double start = mcs2->getTagBitEvt(1<<0, *(dagSteps.begin()+2)+0.1, 0)-82; // dagger first movement ??
        double stop = start + 500;
		double bkgStop = stop + 45;
		
		printf("Start: %f stop: %f diff %f\n", start, stop, stop-start);
		std::vector<input_t> dCts = mcs1->getCoincCounts(
			[start, stop](input_t x)->input_t{x.realtime -= start; return x;},
			[start, stop](input_t x)->bool{return x.realtime > start && x.realtime < stop;}
		);
		std::vector<input_t> bkgdCts = mcs1->getCoincCounts(
			[stop, bkgStop](input_t x)->input_t{return x;},
			[stop, bkgStop](input_t x)->bool{return x.realtime > stop + 5 && x.realtime < bkgStop;}
		);
		if(dCts.size() < 100 || dCts.empty()) {
			return;
		}

		if(dCts.back().realtime > 0 && dCts.back().realtime < 100000) {
			printf("Added - %lu\n", dCts.size());
			printf("Bkg - %lu\n", bkgdCts.size());
			std::for_each(dCts.begin(), dCts.end(), [&dipSumS](input_t x){dipSumS.Fill(x.realtime);});
		}
	};
	
	// Fill our other dip sum histogram (on the Active Cleaner)
	auto acSummer = [&dipSumS](Run* mcs1, Run* mcs2) {
		std::vector<input_t> foilCtsFill = mcs2->getCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return x.ch == 13 && x.realtime < 150.0;}
		);
		if(foilCtsFill.size() < 5000 || foilCtsFill.size() > 60000) {
			return;
		}
		std::vector<input_t> foilCtsHold = mcs2->getCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return x.ch == 13 && x.realtime > 200.0;}
		);
		if(foilCtsHold.size() > 400) {
			return;
		}
		
		std::vector<input_t> acCts = mcs2->getCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return x.ch == 11;}
		);
		std::vector<input_t> coincs = getSelfCoincs(acCts, 50, 5000, 6);
				
		printf("Counts: ");
		std::for_each(coincs.begin(), coincs.end(), [](input_t x){
			if(x.realtime > 210 && x.realtime < 1710) {
				printf("%f,", x.realtime);
			}
		});
		printf("\n");
	};*/

 /* 95(after Run runMCS2(coincWindow, peSumWindow, peSum, runNo,...)):
  			Run run(coincWindow, peSumWindow, peSum, runNo, coincMode, "/media/daq/ssd/2016_2017_data/replayed_data/processed_output_%05d.root"); 
			Run runMCS1(50, 500, 2, runNo, 2, "/media/daq/storage/2016_2017_data/replayed_data/processed_output_%05d.root", "tmcs_0");	
			Run runMCS2(coincWindow, peSumWindow, peSum, runNo, coincMode, "/Volumes/SanDisk/2016-2017/processed_output_%05d.root", "tmcs_1");
			dipSummer(&runMCS1);
			coincACSummer(&runMCS1, &runMCS2);
			ttneAC(&runMCS1, &runMCS2);
			normNByDipSing(&runMCS1);
			bkgRunBkg(&runMCS1);
			noiseFinder(&runMCS1);
			bkgSummer(&runMCS1);
			rayleighPeriodicTest(&runMCS1);
			protheroePeriodicTest(&runMCS1);

/* 99(After that }; at the end of previous code blob):
 		ttne.SaveAs("ttne/ttne.root");
		ttpe.SaveAs("ttne/ttpe.root");
		TH1D dipSumS("summedDipS", "summedDipS", 500, 0, 200);
		TH1D dipSumL("summedDipL", "sumemdDipL", 500, -300, 200);
		std::for_each(sCts.begin(), sCts.end(), [&dipSumS](double x){dipSumS.Fill(x);});
		std::for_each(lCts.begin(), lCts.end(), [&dipSumL](double x){dipSumL.Fill(x);});
		std::for_each(sCts.begin(), sCts.end()-1, [](double x){printf("%f,", x);});
		printf("%f\n", sCts.back());
		std::for_each(lCts.begin(), lCts.end()-1, [](double x){printf("%f,", x);});
		printf("%f\n", lCts.back());
		dipSumS.Sumw2();
		dipSumS.Scale(1.0/dipSumS.Integral());
		dipSumS.SaveAs("summedDaggerTraces/normalizedDipShort.root");
		dipSumS.SaveAs("summedDaggerTraces/summedDipShort.root");
		dipSumL.SaveAs("summedDaggerTraces/summedDipLong.root");
		singleBkg.SaveAs("summedDaggerTraces/summedLongCtsSingle.root");*/
/* 101(before the last Return 0):
	hand.foreach(fitFill);
	hand.foreach(normNByDip);
	saveAllPMTHits(pmtAHits);
	
	hand.foreach(writeCoincHist);
	summedHist.SaveAs("summedDaggerTraces/summedTrace.root");
	pmtAWaveform.SaveAs("waveforms/waveformA.root");
	pmtBWaveform.SaveAs("waveforms/waveformB.root");
	printf("Total # coinc: %ld\n", totalNumCoinc);
	sumphsA.SaveAs("pulseHeightSpectra/phsA.root");
	sumphsB.SaveAs("pulseHeightSpectra/phsB.root");
*/
	/*
			// Fit the summed data
			auto maxhgx = max_element(std::begin(numHGXhits),std::end(numHGXhits));
			auto minhgx = min_element(std::begin(numHGXhits),std::end(numHGXhits));
			printf("%d\n",maxhgx);
			if (maxhgx == minhgx) {
				printf("Fitting summed fills !!! ");
				printf("Using hardcoded H-GX spacings (5.4s) !!!\n");
				ExpFillFree func;
				for(int it = 0; it < 55; it++) {
					/*if(*it+3.0 > endFillTime) {
						beamHits.pop_back();
						continue;
					}*/
				/*	func.addOffset((double)it*5.425);
				}
				printf("%d",func.getNumOffsets() );
				// Define ROOT fitting tools
				//TF1* fitd1 = new TF1("fitd1", func, 0.0, endFillTime, 4);
				//TF1* fitd2 = new TF1("fitd2", func, 0.0, endFillTime, 4);
				//TF1* fitol = new TF1("fitol", func, 0.0, endFillTime, 4);
				TF1* fitba = new TF1("fitba", func, 0.0, endFillTime, 4);
				//TF1* fitsp = new TF1("fitsp", func, 0.0, endFillTime, 4);
				
				// Set a stupid amount of parameters
				fitba->SetParameter(0, 3.0);
				fitba->SetParLimits(0, 0.0, 5.00);
				fitba->SetParameter(1, 30000);
				fitba->SetParLimits(1, 0.0, 99999);
				fitba->SetParameter(2, 0.5);
				//fitba->SetParLimits(2, 0.0, 100.0);
				fitba->SetParLimits(2, 0.49, 0.51);
				fitba->SetParameter(3, 20.0);
				fitba->SetParLimits(3, 19.9, 20.1);
				fillSumBa.Fit(fitba);
				
				//for (int det = 1; det < 6; det++) {
					printf("FillData - %d,%d,%f,%f\n", 4, -1, fitba->GetChisquare(), (double)fitba->GetNDF());
					for (int par = 0; par < 4; par++) {
						printf("FillData - %d,%d,%f,%f\n", 4, par, fitba->GetParameter(par), fitba->GetParError(par));
					}
				//}
				
			} else {
				printf("Not fitting summed fills !!!");
				printf(" Inconsistent number of H-GX pulses !!!\n");
			}
			
			// Save histograms
			
			*///} break;
		
//		default : 
//			printf("Not saving external ROOT files!\n");
//	}
	
