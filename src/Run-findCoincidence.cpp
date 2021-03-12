#include "../inc/Run.hpp"
#define NANOSECOND .000000001

/*	------------------------------------------------------------------------------------------------
	Author: Nathan B. Callahan
	Editor: Frank M. Gonzalez
	
	This function takes the data and (empty) coincidence vectors, and the coincidence cuts and fills
	the coincidence vector with the coincidences. It also makes a few auxiliarry plots of the
	spectrum and relative timing of coincidences.
	
	The function iterates through all entries of the data vector. It begins coincidence finding by
	searching for an event in the betas. If it finds a beta event, then it begins a local search
	forward in time bound by the upper timing window for a coincident gamma.
	
	If the current beta event was within a previous beta event's timing window, then reject the
	event and continue.
	
	If a coincident beta is found, the search is terminated and the event is rejected.
	Any gamma of w/ energy > 250 keV is added to a count of gammas in the window.
	If a gamma is inside of the window and inside of the energy cuts, then the number of coincident
	gammas is incremented. The auxiliary plots are filled at this time.
	
	If there is only one gamma in the window and it met the timing and energy cuts, then the event
	is placed into the coincidence vector.
	------------------------------------------------------------------------------------------------	*/

void Run::findCoincidenceFixedTele() {
	
	// check the data table and load it into ROOT tree 
	if(data.empty()) {
		return;
	}

	// initialize iterators and variables 
	int i;
	int cur = 0;
	int tailIt = 0;
	//int waveIt = 0;
	
	int ch1 = coinChan1 + mcsOffset;
	int ch2 = coinChan2 + mcsOffset;
	
	// initialize variables for our pmt hits
	int ch1PESum;
	int ch1InstSum;
	int ch2PESum;
	int ch2InstSum;
	int promptSum;
	std::vector<input_t> pmtAHits;
	std::vector<input_t> pmtBHits;
	std::vector<input_t> coincHits;
	std::vector<long> coincIndices;
	coinc_t cParam;
	cParam.rate = 0; // Rate of each event is 0
	double tailT;
	
	unsigned long nCoinc = 0;
	// look at each entry of our data vector.
	for(i = 0; i < data.size(); i++) {
		
		// Clear previous buffers 
		ch1PESum = 0;
		ch2PESum = 0;
		promptSum = 0;
		pmtAHits.clear();
		pmtBHits.clear();
		coincHits.clear();
		tailT = 0.0;
		input_t prevEvt;
		
		// Make sure element i is in the dagger
		if(data.at(i).ch != ch1 && data.at(i).ch != ch2) { continue; }
		
		// Load info into our data hit vectors
		if(data.at(i).ch == ch1) {
			ch1PESum += 1;
			pmtAHits.push_back(data.at(i));
			coincHits.push_back(data.at(i));
		}
		if(data.at(i).ch == ch2) {
			ch2PESum += 1;
			pmtBHits.push_back(data.at(i));
			coincHits.push_back(data.at(i));
		}
		// Once we load our dataset, search forwards to find coincidences
		for(cur = i+1; cur < data.size(); cur++) { 
			
			// Check the times of our two paired events. 
			// If the times are not about the same, break. 
			// We don't have a coincidence! 
			if(data.at(cur).time - data.at(i).time > (unsigned long) (coincWindow/0.8)) {
				break; 
			}
			// Make sure element cur is also in the dagger
			if(data.at(cur).ch != ch1 && data.at(cur).ch != ch2) { continue; }
			
			// Count our coincidences in the same channel i (self Coinc)
			if(data.at(cur).ch == data.at(i).ch) {
				if(data.at(cur).ch == ch1) {
					ch1PESum += 1;
					promptSum += 1;
					pmtAHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == ch2) {
					ch2PESum += 1;
					promptSum += 1;
					pmtBHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
			}
			
			// Count our coincidences in different channels
			if(data.at(cur).ch != data.at(i).ch) {				
				if(data.at(cur).ch == ch1) {
					ch1PESum += 1;
					promptSum += 1;
					pmtAHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == ch2) {
					ch2PESum += 1;
					promptSum += 1;
					pmtBHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
				prevEvt = data.at(cur);
				// integrate the tail end. Add data to the sums of the 
				// two channels. 
				for(tailIt = cur+1; tailIt < data.size(); tailIt++) {
					//if(data.at(tailIt).time - data.at(i).time > (unsigned long) (peSumWindow/0.8)) {
					if(data.at(tailIt).time - prevEvt.time > (unsigned long) (peSumWindow/0.8)) {
						break;
					}
					if(data.at(tailIt).ch == ch1) {
						ch1PESum += 1;
						pmtAHits.push_back(data.at(tailIt));
						coincHits.push_back(data.at(tailIt));
						// Add the prompt sum for the first half of the tail
						if((data.at(tailIt).time - data.at(i).time) < (unsigned long) (promptWindow/0.8)) {
							promptSum += 1;
						}
					}
					if(data.at(tailIt).ch == ch2) {
						ch2PESum += 1;
						pmtBHits.push_back(data.at(tailIt));
						coincHits.push_back(data.at(tailIt));
						if((data.at(tailIt).time - data.at(i).time) < (unsigned long) (promptWindow/0.8)) {
							promptSum += 1;
						}
					}
					
					//tailT = data.at(tailIt).realtime; // Count the last time of tail
					prevEvt = data.at(tailIt);
				}
				// Check to see if we found a neutron
				// This particular one puts a 1.5*deadtime 
				if((promptSum >= peSum) &&
					(prevEvt.realtime - data.at(i).realtime > (double)(this->fastDT*(ch1PESum+ch2PESum)*0.8*NANOSECOND))) {
					// Add to c++ vector
					// Set cParam here.
					cParam.realtime = data.at(i).realtime;
					cParam.time = data.at(i).time;
					cParam.pmt1 = ch1PESum;
					cParam.pmt2 = ch2PESum;
					cParam.prompt = promptSum;
					cParam.length = prevEvt.realtime - data.at(i).realtime + peSumWindow*NANOSECOND;
					//cParam.length = peSumWindow*NANOSECOND;
					cParam.index = nCoinc;
					nCoinc += 1;
					// Calculate mean times for pmtA and pmtB
					double meanA = 0.0;
					if (pmtAHits.size() > 1) {
						for (size_t jj = 1; jj < pmtAHits.size(); jj++) {
							meanA += pmtAHits.at(jj).realtime - pmtAHits.at(jj-1).realtime;
						}
					}
					cParam.avg1 = meanA / (double)pmtAHits.size();
					double meanB = 0.0;
					if (pmtBHits.size() > 1) {
						for (size_t jj = 1; jj < pmtBHits.size(); jj++) {
							meanB += pmtBHits.at(jj).realtime - pmtBHits.at(jj-1).realtime;
						}
					}
					cParam.avg2 = meanB / (double)pmtBHits.size();				
					
					//coinc.push_back(data.at(i));
					coinc.push_back(cParam);
					coincIndices.push_back(i);
					pmtACoincHits.push_back(pmtAHits);
					pmtBCoincHits.push_back(pmtBHits);
					for (size_t longTail = tailIt; longTail < data.size(); longTail++) {
						
						if (data.at(longTail).realtime - data.at(i).realtime > 0.001) {
							break;
						}
						
						if (data.at(longTail).ch == ch1 || data.at(longTail).ch == ch2) {
							coincHits.push_back(data.at(longTail));
						}
						//printf("%d,%d,%d%,u\n",i,cur,tailIt,longTail);
					}
					//}
					pmtCCoincHits.push_back(coincHits);
					i = tailIt-1; 
					cur = tailIt-1;
				
					break;
				}
				break;
			}
		}
	}
}

/* Coincidence timer -- subset of run. 
 * Fixed means we don't move the timing window. */
void Run::findCoincidenceFixed() {
	
	// check the data table and load it into ROOT tree 
	if(data.empty()) {
		return;
	}

	// initialize iterators and variables 
	int i;
	int cur = 0;
	int tailIt = 0;
	//int waveIt = 0;
	
	int ch1 = coinChan1 + mcsOffset;
	int ch2 = coinChan2 + mcsOffset;
	
	// initialize variables for our pmt hits
	int ch1PESum;
	int ch2PESum;
	int promptSum;
	std::vector<input_t> pmtAHits;
	std::vector<input_t> pmtBHits;
	std::vector<input_t> coincHits;
	std::vector<long> coincIndices;
	coinc_t cParam;
	cParam.rate = 0; // Rate of each event is 0
	double tailT;
	
	unsigned long nCoinc = 0;
	// look at each entry of our data vector.
	for(i = 0; i < data.size(); i++) {
		
		// Clear previous buffers 
		ch1PESum = 0;
		ch2PESum = 0;
		promptSum = 0;
		pmtAHits.clear();
		pmtBHits.clear();
		coincHits.clear();
		tailT = 0.0;

		// Make sure element i is in the dagger
		if(data.at(i).ch != ch1 && data.at(i).ch != ch2) { continue; }
		
		// Load info into our data hit vectors
		if(data.at(i).ch == ch1) {
			ch1PESum += 1;
			pmtAHits.push_back(data.at(i));
			coincHits.push_back(data.at(i));
		}
		if(data.at(i).ch == ch2) {
			ch2PESum += 1;
			pmtBHits.push_back(data.at(i));
			coincHits.push_back(data.at(i));
		}
		
		// Once we load our dataset, search forwards to find coincidences
		for(cur = i+1; cur < data.size(); cur++) { 
			
			// Check the times of our two paired events. 
			// If the times are not about the same, break. 
			// We don't have a coincidence! 
			if(data.at(cur).time - data.at(i).time > (unsigned long) (coincWindow/0.8)) {
				break; 
			}
			// Make sure element cur is also in the dagger
			if(data.at(cur).ch != ch1 && data.at(cur).ch != ch2) { continue; }
			
			// Count our coincidences in the same channel i (self Coinc)
			if(data.at(cur).ch == data.at(i).ch) {
				if(data.at(cur).ch == ch1) {
					ch1PESum += 1;
					promptSum += 1;
					pmtAHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == ch2) {
					ch2PESum += 1;
					promptSum += 1;
					pmtBHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
			}
			
			// Count our coincidences in different channels
			if(data.at(cur).ch != data.at(i).ch) {				
				if(data.at(cur).ch == ch1) {
					ch1PESum += 1;
					promptSum += 1;
					pmtAHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == ch2) {
					ch2PESum += 1;
					promptSum += 1;
					pmtBHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
				// integrate the tail end. Add data to the sums of the 
				// two channels. 
				for(tailIt = cur+1; tailIt < data.size(); tailIt++) {
					if(data.at(tailIt).time - data.at(i).time > (unsigned long) (peSumWindow/0.8)) {
						break;
					}
					if(data.at(tailIt).ch == ch1) {
						ch1PESum += 1;
						pmtAHits.push_back(data.at(tailIt));
						coincHits.push_back(data.at(tailIt));
						// Add the prompt sum for the first half of the tail
						if((data.at(tailIt).realtime - data.at(i).realtime) < peSumWindow/2.) {
							promptSum += 1;
						}
					}
					if(data.at(tailIt).ch == ch2) {
						ch2PESum += 1;
						pmtBHits.push_back(data.at(tailIt));
						coincHits.push_back(data.at(tailIt));
						if((data.at(tailIt).realtime - data.at(i).realtime) < peSumWindow/2.) {
							promptSum += 1;
						}
					}
					
					tailT = data.at(tailIt).realtime; // Count the last time of tail
				}
				// Check to see if we found a neutron
				if(ch1PESum + ch2PESum >= peSum) {
					// Add to c++ vector
					// Set cParam here.
					cParam.realtime = data.at(i).realtime;
					cParam.time = data.at(i).time;
					cParam.pmt1 = ch1PESum;
					cParam.pmt2 = ch2PESum;
					cParam.prompt = promptSum;
					cParam.length = peSumWindow*NANOSECOND;
					cParam.index = nCoinc;
					nCoinc += 1;
					
					// Calculate mean times for pmtA and pmtB
					double meanA = 0.0;
					if (pmtAHits.size() > 1) {
						for (size_t jj = 1; jj < pmtAHits.size(); jj++) {
							meanA += pmtAHits.at(jj).realtime - pmtAHits.at(jj-1).realtime;
						}
					}
					cParam.avg1 = meanA / (double)pmtAHits.size();
					double meanB = 0.0;
					if (pmtBHits.size() > 1) {
						for (size_t jj = 1; jj < pmtBHits.size(); jj++) {
							meanB += pmtBHits.at(jj).realtime - pmtBHits.at(jj-1).realtime;
						}
					}
					cParam.avg2 = meanB / (double)pmtBHits.size();				
					
					//coinc.push_back(data.at(i));
					coinc.push_back(cParam);
					coincIndices.push_back(i);
					pmtACoincHits.push_back(pmtAHits);
					pmtBCoincHits.push_back(pmtBHits);
					pmtCCoincHits.push_back(coincHits);	
					i = tailIt-1; 
					cur = tailIt-1;
					break;
				}
				break;
			}
		}
	}
}

/* Another coincidence timer. The difference between this one and the 
 * other type is that this one keeps track of the previous event. */
void Run::findCoincidenceMoving() {
	
	// check the data table and load it into the ROOT tree
	/*if(data.empty()) {
		this->readDataRoot();
	}*/
	if(data.empty()) {
		return;
	}

	// initialize iterators and variables
	int i;
	int cur = 0;
	int tailIt = 0;
	//int waveIt = 0;
	
	int ch1 = coinChan1 + mcsOffset;
	int ch2 = coinChan2 + mcsOffset;
	// initialize variables for our PMT hits
	int ch1PESum;
	int ch2PESum;
	int promptSum;
	std::vector<input_t> pmtAHits;
	std::vector<input_t> pmtBHits;
	std::vector<input_t> coincHits;
	std::vector<long> coincIndices;
	coinc_t cParam;
	cParam.rate = 0; // Rate of each event is zero unless you do something fancy later
	
	unsigned long nCoinc = 0;
	// look at each entry of the data vector
	for(i = 0; i < data.size(); i++) {
		
		// clear previous buffers 
		ch1PESum = 0;
		ch2PESum = 0;
		promptSum = 0;
		pmtAHits.clear();
		pmtBHits.clear();
		coincHits.clear();
		input_t prevEvt;
		
		// for coincidence measurements we only want dagger hits
		if(data.at(i).ch != ch1 && data.at(i).ch != ch2) { continue; }
		
		// load info into our channel hit vectors
		if(data.at(i).ch == ch1) {
			ch1PESum += 1;
			promptSum += 1;
			pmtAHits.push_back(data.at(i));
			coincHits.push_back(data.at(i));
		}
		if(data.at(i).ch == ch2) {
			ch2PESum += 1;
			promptSum += 1;
			pmtBHits.push_back(data.at(i));
			coincHits.push_back(data.at(i));
		}
		
		// once we've loaded our dataset, search forwards to find coincidence
		for(cur = i+1; cur < data.size(); cur++) {
			
			// check the times of our two events. If the two are too far
			// apart, break because it's not a coincidence!
			//unsigned long ulcast = 1/(0.8*NANOSECOND);
			if(data.at(cur).time - data.at(i).time > (unsigned long)(coincWindow/0.8)) {
			//if(data.at(cur).realtime - data.at(i).realtime > coincWindow*NANOSECOND) {
				break;
			}
			// only count coincidences in the dagger
			if(data.at(cur).ch != ch1 && data.at(cur).ch != ch2) { continue; }
			// count coincidences if we find a second event in channel i
			if(data.at(cur).ch == data.at(i).ch) {
				if(data.at(i).ch == ch1) {
					ch1PESum  += 1;
					promptSum += 1;
					pmtAHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
				if(data.at(i).ch == ch2) {
					ch2PESum  += 1;
					promptSum += 1;
					pmtBHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
			}
			// count our coincidences on different channels
			if(data.at(cur).ch != data.at(i).ch) {
				if(data.at(cur).ch == ch1) {
					ch1PESum  += 1;
					promptSum += 1;
					pmtAHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == ch2) {
					ch2PESum  += 1;
					promptSum += 1;
					pmtBHits.push_back(data.at(cur));
					coincHits.push_back(data.at(cur));
				}
				
				// save the previous event as a new data point
				prevEvt = data.at(cur);
				
				// Integrate the tail end.
				// Add counts into the right channel 
				for(tailIt = cur+1; tailIt < data.size(); tailIt++) {
					if((data.at(tailIt).time - prevEvt.time) > (unsigned long)(peSumWindow/0.8)) {
					//if((data.at(tailIt).realtime - prevEvt.realtime) > peSumWindow*NANOSECOND) {
						break;
					}
					if(data.at(tailIt).ch != ch1 && data.at(tailIt).ch != ch2) { continue; }
					if(data.at(tailIt).ch == ch1) {
						ch1PESum += 1;
						pmtAHits.push_back(data.at(tailIt));
						coincHits.push_back(data.at(tailIt));
						if((data.at(tailIt).time - data.at(i).time) < (unsigned long)(promptWindow/0.8)) {
							promptSum += 1;
						}
					}
					if(data.at(tailIt).ch == ch2) {
						ch2PESum += 1;
						pmtBHits.push_back(data.at(tailIt));
						coincHits.push_back(data.at(tailIt));
						if((data.at(tailIt).time - data.at(i).time) < (unsigned long)(promptWindow/0.8)) {
							promptSum += 1;
						}
					}
					prevEvt = data.at(tailIt);
				}
				// check to see if we've found a neutron!
				if(promptSum >= peSum) {
				//if(ch1PESum + ch2PESum >= peSum) {
					// fill our vectors
					// Set cParam here.
					cParam.realtime = data.at(i).realtime;
					cParam.time = data.at(i).time;
					cParam.pmt1 = ch1PESum;
					cParam.pmt2 = ch2PESum;
					cParam.prompt = promptSum;
					cParam.length = prevEvt.realtime - data.at(i).realtime;
					cParam.index = nCoinc;
					nCoinc += 1;
					// Calculate mean times for pmtA and pmtB
					double meanA = 0.0;
					if (pmtAHits.size() > 1) {
						for (size_t jj = 1; jj < pmtAHits.size(); jj++) {
							meanA += pmtAHits.at(jj).realtime - pmtAHits.at(jj-1).realtime;
						}
						cParam.avg1 = meanA / ((double)pmtAHits.size() - 1.0);
					} else {
						cParam.avg1 = 0.0;
					}
					
					double meanB = 0.0;
					if (pmtBHits.size() > 1) {
						for (size_t jj = 1; jj < pmtBHits.size(); jj++) {
							meanB += pmtBHits.at(jj).realtime - pmtBHits.at(jj-1).realtime;
						}
						cParam.avg2 = meanB / ((double)pmtBHits.size() - 1.0);
					} else {
						cParam.avg2 = 0.0;
					}
										
					coinc.push_back(cParam);
					//coinc.push_back(data.at(i));
					coincIndices.push_back(i);
					pmtACoincHits.push_back(pmtAHits);
					pmtBCoincHits.push_back(pmtBHits);
					for (int longTail = tailIt; longTail < data.size(); longTail++) {
						if (data.at(longTail).realtime - data.at(i).realtime > 0.001) {
							break;
						}
						
						if (data.at(longTail).ch == ch1 || data.at(longTail).ch == ch2) {
							coincHits.push_back(data.at(longTail));
						}
					}
					pmtCCoincHits.push_back(coincHits);
					//phsA.Fill(ch1PESum);
					//phsB.Fill(ch2PESum);
					
					i = tailIt-1; 
					cur = tailIt-1;
					break;
				}
				break;
			}
		}
	}
}

/* Another coincidence timer. This one doesn't require multiple PMTs
 * to be hit in the initial window. */
void Run::findCoincidenceNoPMT() {
	
	// check the data table and load it into the ROOT tree
	if(data.empty()) {
		return;
	}

	// initialize iterators and variables
	int i;
	int cur = 0;
	int tailIt = 0;
	//int waveIt = 0;
	
	int ch1 = coinChan1 + mcsOffset;
	int ch2 = coinChan2 + mcsOffset;
	// initialize variables for our PMT hits
	int ch1PESum;
	int ch2PESum;
	int promptSum;
	std::vector<input_t> pmtAHits;
	std::vector<input_t> pmtBHits;
	std::vector<long> coincIndices;
	coinc_t cParam;
	cParam.rate = 0; // Rate of each event starts at zero
	// look at each entry of the data vector
	for(i = 0; i < data.size(); i++) {
		
		// clear previous buffers 
		ch1PESum = 0;
		ch2PESum = 0;
		pmtAHits.clear();
		pmtBHits.clear();
		input_t prevEvt;
		
		// for coincidence measurements we only want dagger hits
		if(data.at(i).ch != ch1 && data.at(i).ch != ch2) { continue; }
		
		// load info into our channel hit vectors
		if(data.at(i).ch == ch1) {
			ch1PESum  += 1;
			promptSum += 1;
			pmtAHits.push_back(data.at(i));
		}
		if(data.at(i).ch == ch2) {
			ch2PESum  += 1;
			promptSum += 1;
			pmtBHits.push_back(data.at(i));
		}
		
		// once we've loaded our dataset, search forwards to find coincidence
		for(cur = i+1; cur < data.size(); cur++) {
			
			// check the times of our two events. If the two are too far
			// apart, break because it's not a coincidence!
			if(data.at(cur).time - data.at(i).time > (unsigned long) (coincWindow/0.8)) {
				break;
			}
			// only count coincidences in the dagger
			if(data.at(cur).ch != ch1 && data.at(cur).ch != ch2) { continue; }
						
			// count our coincidences on different channels
			if(data.at(cur).ch == ch1) {
				ch1PESum  += 1;
				promptSum += 1;
				pmtAHits.push_back(data.at(cur));
			}
			if(data.at(cur).ch == ch2) {
				ch2PESum  += 1;
				promptSum += 1;
				pmtBHits.push_back(data.at(cur));
			}
			
			// save the previous event as a new data point
			prevEvt = data.at(cur);
			
			// Integrate the tail end.
			// Add counts into the right channel 
			for(tailIt = cur+1; tailIt < data.size(); tailIt++) {
				if((data.at(tailIt).time - prevEvt.time) > (unsigned long) (peSumWindow/0.8)) {
					break;
				}
				if(data.at(tailIt).ch != ch1 && data.at(tailIt).ch != ch2) { continue; }
				if(data.at(tailIt).ch == ch1) {
					ch1PESum += 1;
					pmtAHits.push_back(data.at(tailIt));
					if((data.at(tailIt).realtime - data.at(i).realtime) < peSumWindow) {
						promptSum += 1;
					}
				}
				if(data.at(tailIt).ch == ch2) {
					ch2PESum += 1;
					pmtBHits.push_back(data.at(tailIt));
					if((data.at(tailIt).realtime - data.at(i).realtime) < peSumWindow) {
						promptSum += 1;
					}
				}
				prevEvt = data.at(tailIt);
			}
			// check to see if we've found a neutron!
			if(ch1PESum + ch2PESum >= peSum) {
				// fill our vectors
				// Set cParam here.
				cParam.realtime = data.at(i).realtime;
				cParam.time = data.at(i).time;
				cParam.pmt1 = ch1PESum;
				cParam.pmt2 = ch2PESum;
				cParam.prompt = promptSum;
				cParam.length = prevEvt.realtime - data.at(i).realtime;
				
				// Calculate mean times for pmtA and pmtB
				double meanA = 0.0;
				if (pmtAHits.size() > 1) {
					for (size_t jj = 1; jj < pmtAHits.size(); jj++) {
						meanA += pmtAHits.at(jj).realtime - pmtAHits.at(jj-1).realtime;
					}
					cParam.avg1 = meanA / ((double)pmtAHits.size() - 1.0);
				} else {
					cParam.avg1 = 0.0;
				}
				
				double meanB = 0.0;
				if (pmtBHits.size() > 1) {
					for (size_t jj = 1; jj < pmtBHits.size(); jj++) {
						meanB += pmtBHits.at(jj).realtime - pmtBHits.at(jj-1).realtime;
					}
					cParam.avg2 = meanB / ((double)pmtBHits.size() - 1.0);
				} else {
					cParam.avg2 = 0.0;
				}
									
				coinc.push_back(cParam);
				//coinc.push_back(data.at(i));
				coincIndices.push_back(i);
				pmtACoincHits.push_back(pmtAHits);
				pmtBCoincHits.push_back(pmtBHits);
				//phsA.Fill(ch1PESum);
				//phsB.Fill(ch2PESum);
				
				i = tailIt-1; 
				cur = tailIt-1;
				break;
			}
		}
	}
}


/* Active cleaner coincidence timer
 * This is coincidence type 1 for AC
 */
void Run::findCoincidenceAC() {
	
	// check the data table and load it into the ROOT tree
	/*if(data.empty()) {
		this->readDataRoot();
	}*/
	if(data.empty()) {
		return;
	}

	// initialize iterators and variables
	int i;
	int cur = 0;
	int tailIt = 0;
	//int waveIt = 0;
	
	// initialize variables for our PMT hits
	int ch1PESum;
	int ch2PESum;
	int promptSum;
	int acPESum;
	std::vector<input_t> pmtAHits;
	std::vector<input_t> pmtBHits;
	std::vector<input_t> acHits;
	std::vector<long> coincIndices;
	coinc_t cParam;
	double tailT;
	
	// look at each entry of the data vector
	for(i = 0; i < data.size(); i++) {
		
		// clear previous buffers 
		ch1PESum = 0;
		ch2PESum = 0;
		promptSum = 0;
		acPESum  = 0;
		pmtAHits.clear();
		pmtBHits.clear();
		acHits.clear();
		
		input_t prevEvt;
		
		// for coincidence measurements we only want dagger hits
		if(data.at(i).ch != 11 && (data.at(i).ch != 14 && data.at(i).ch != 15)) { continue; }
		
		// load info into our channel hit vectors
		if(data.at(i).ch == 11) {
			acPESum += 1;
			acHits.push_back(data.at(i));
		}
		if(data.at(i).ch == 14) {
			ch1PESum  += 1;
			promptSum += 1;
			pmtAHits.push_back(data.at(i));
		}
		if(data.at(i).ch == 15) {
			ch2PESum  += 1;
			promptSum += 1;
			pmtBHits.push_back(data.at(i));
		}
		
		// once we've loaded our dataset, search forwards to find coincidence
		for(cur = i+1; cur < data.size(); cur++) {
			
			// check the times of our two events. If the two are too far
			// apart, break because it's not a coincidence!
			if(data.at(cur).time - data.at(i).time > (unsigned long) (coincWindow/0.8)) {
				break;
			}
			// only count coincidences in the dagger or AC
			if(data.at(cur).ch != 11 && (data.at(cur).ch != 14 && data.at(cur).ch != 15)) { continue; }
			// count coincidences if we find a second event in channel i
			if(data.at(cur).ch == data.at(i).ch) {
				if(data.at(i).ch == 11) {
					acPESum += 1;
					acHits.push_back(data.at(cur));
				}
				if(data.at(i).ch == 14) {
					ch1PESum  += 1;
					promptSum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(i).ch == 15) {
					ch2PESum  += 1;
					promptSum += 1;
					pmtBHits.push_back(data.at(cur));
				}
			}
			// count our coincidences on different channels
			if(data.at(cur).ch != data.at(i).ch) {
				if(data.at(cur).ch == 11) {
					acPESum += 1;
					acHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == 14) {
					ch1PESum  += 1;
					promptSum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == 15) {
					ch2PESum  += 1;
					promptSum += 1;
					pmtBHits.push_back(data.at(cur));
				}
				
				// Force AC to trigger
				if (acPESum == 0) {
					continue;
				}
				
				// save the previous event as a new data point
				prevEvt = data.at(cur);
				
				// integrate the tail end. Add counts into the right 
				// channel 
				for(tailIt = cur+1; tailIt < data.size(); tailIt++) {
					if((data.at(tailIt).time - prevEvt.time) > (unsigned long) (peSumWindow/0.8)) {
						break;
					}
					if(data.at(tailIt).ch == 11 && (data.at(tailIt).ch != 14 && data.at(tailIt).ch != 15)) { continue; }
					if(data.at(tailIt).ch == 11) {
						acPESum += 1;
						acHits.push_back(data.at(tailIt));
					}
					if(data.at(tailIt).ch == 14) {
						ch1PESum += 1;
						if((data.at(tailIt).realtime - data.at(i).realtime) < peSumWindow) {
							promptSum += 1;
						}
						pmtAHits.push_back(data.at(tailIt));
					}
					if(data.at(tailIt).ch == 15) {
						ch2PESum += 1;
						if((data.at(tailIt).realtime - data.at(i).realtime) < peSumWindow) {
							promptSum += 1;
						}
						pmtBHits.push_back(data.at(tailIt));
					}
					prevEvt = data.at(tailIt);
				}
				
				if ((ch1PESum + acPESum >= peSum) || (ch2PESum + acPESum >= peSum)) {
					// fill our vectors
					// Set cParam here.
					cParam.realtime = data.at(i).realtime;
					cParam.time = data.at(i).time;
					cParam.pmt1 = ch1PESum;
					cParam.pmt2 = ch2PESum;
					cParam.prompt = acPESum; // sliding AC counts in here.
					cParam.length = prevEvt.realtime - data.at(i).realtime;
					
					coinc.push_back(cParam);
					//coinc.push_back(data.at(i));
					coincIndices.push_back(i);
					pmtACoincHits.push_back(pmtAHits);
					pmtBCoincHits.push_back(pmtBHits);
					acCoincHits.push_back(acHits);
					//phsA.Fill(ch1PESum);
					//phsB.Fill(ch2PESum);
					// put deadtime on counted neutrons
					i = tailIt-1; 
					cur = tailIt-1;
					break;
				}
				break;
			}
		}
	}
}

/* AC Anti-coincidence algorithm
 * Looks for coincidences and only returns 
 * ANTI-coincidences, i.e. AC hits with no coincidence 
 * going on alongside them. */
void Run::findAntiCoincidenceAC() {
	
	// check the data table and load it into the ROOT tree
	/*if(data.empty()) {
		this->readDataRoot();
	}*/
	if(data.empty()) {
		return;
	}

	// initialize iterators and variables
	int i;
	int cur = 0;
	int tailIt = 0;
	//int waveIt = 0;
	
	// initialize variables for our PMT hits
	int ch1PESum;
	int ch2PESum;
	int acPESum;
	std::vector<input_t> pmtAHits;
	std::vector<input_t> pmtBHits;
	std::vector<input_t> acHits;
	std::vector<long> coincIndices;

	coinc_t cParam;
	// look at each entry of the data vector
	for(i = 0; i < data.size(); i++) {
		
		// clear previous buffers 
		ch1PESum = 0;
		ch2PESum = 0;
		acPESum  = 0;
		pmtAHits.clear();
		pmtBHits.clear();
		acHits.clear();
		
		input_t prevEvt;
		
		// for coincidence measurements we only want dagger hits
		if(data.at(i).ch != 11 && (data.at(i).ch != 14 && data.at(i).ch != 15)) { continue; }
		
		// load info into our channel hit vectors
		if(data.at(i).ch == 11) {
			acPESum += 1;
			acHits.push_back(data.at(i));
		}
		if(data.at(i).ch == 14) {
			ch1PESum += 1;
			pmtAHits.push_back(data.at(i));
		}
		if(data.at(i).ch == 15) {
			ch2PESum += 1;
			pmtBHits.push_back(data.at(i));
		}
		
		// once we've loaded our dataset, search forwards to find coincidence
		for(cur = i+1; cur < data.size(); cur++) {
			
			// check the times of our two events. If the two are too far
			// apart, break because it's not a coincidence!
			if(data.at(cur).time - data.at(i).time > (unsigned long) (coincWindow/0.8)) {
				// If data.at(i) is an AC hit, we should add it!
				if ((data.at(i).ch == 11) && (data.at(cur).ch == 11)) {
					
					// Set cParam here.
					cParam.realtime = data.at(i).realtime;
					cParam.time = data.at(i).time;
					cParam.pmt1 = 0;
					cParam.pmt2 = 0;
					cParam.length = data.at(cur).realtime - data.at(i).realtime;
					
					coinc.push_back(cParam);
					
					//coinc.push_back(data.at(i));
					coincIndices.push_back(i);
					acCoincHits.push_back(acHits);
				}
				break;
			}
			// only count coincidences in the dagger or AC
			if(data.at(cur).ch != 11 && (data.at(cur).ch != 14 && data.at(cur).ch != 15)) { continue; }
			//if(data.at(cur).ch != 14 && data.at(cur).ch != 15) { continue; }
			// count coincidences if we find a second event in channel i
			if(data.at(cur).ch == data.at(i).ch) {
				if(data.at(i).ch == 11) {
					acPESum += 1;
					acHits.push_back(data.at(cur));
				}
				if(data.at(i).ch == 14) {
					ch1PESum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(i).ch == 15) {
					ch2PESum += 1;
					pmtBHits.push_back(data.at(cur));
				}
			}
			// count our coincidences on different channels
			if(data.at(cur).ch != data.at(i).ch) {
				if(data.at(cur).ch == 11) {
					acPESum += 1;
					acHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == 14) {
					ch1PESum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == 15) {
					ch2PESum += 1;
					pmtBHits.push_back(data.at(cur));
				}
				
				// save the previous event as a new data point
				prevEvt = data.at(cur);
				
				// integrate the tail end. Add counts into the right 
				// channel 
				for(tailIt = cur+1; tailIt < data.size(); tailIt++) {
					if((data.at(tailIt).time - prevEvt.time) > (unsigned long) (peSumWindow/0.8)) {
						break;
					}
					if(data.at(tailIt).ch != 11 && (data.at(tailIt).ch != 14 && data.at(tailIt).ch != 15)) { continue; }
					//if(data.at(tailIt).ch != 14 && data.at(tailIt).ch != 15) { continue; }
					if(data.at(tailIt).ch == 11) {
						acPESum += 1;
						acHits.push_back(data.at(tailIt));
					}
					if(data.at(tailIt).ch == 14) {
						ch1PESum += 1;
						pmtAHits.push_back(data.at(tailIt));
					}
					if(data.at(tailIt).ch == 15) {
						ch2PESum += 1;
						pmtBHits.push_back(data.at(tailIt));
					}
					prevEvt = data.at(tailIt);
				}
				// check to see if we've found a neutron!
				if (ch1PESum + ch2PESum + acPESum >= peSum) {
					i = tailIt-1; 
					cur = tailIt-1;
					break;
					// fill our vectors
					/*coinc.push_back(data.at(i));
					coincIndices.push_back(i);
					acCoincHits.push_back(acHits);*/
					//phsA.Fill(ch1PESum);
					//phsB.Fill(ch2PESum);
					
					
				}
				break;
			}
		}
	}
}

/* Coincidence timer -- subset of run. 
 * Fixed means we don't move the timing window. 
 * This generates ROOT histograms */
/*void Run::findFixedHistos() {
	
	// check the data table and load it into ROOT tree 
	if(data.empty()) {
		this->readDataRoot();
	}
	if(data.empty()) {
		return;
	}
	
	// Load coincidence data
	// initialize iterators and variables 
	int i;
	int cur = 0;
	int tailIt = 0;
	int waveIt = 0;
	
	// initialize variables for our pmt hits
	int ch1PESum;
	int ch2PESum;
	std::vector<input_t> pmtAHits;
	std::vector<input_t> pmtBHits;
	std::vector<long> coincIndices;

	// look at each entry of our data vector.
	for(i = 0; i < data.size(); i++) {
		
		// clear previous buffers 
		ch1PESum = 0;
		ch2PESum = 0;
		pmtAHits.clear();
		pmtBHits.clear();

		// for coincidence measurement, we only want dagger hits
		if(data.at(i).ch != 1 && data.at(i).ch != 2) { continue; }
		
		// load info into our data hit vectors
		if(data.at(i).ch == 1) {
			ch1PESum += 1;
			pmtAHits.push_back(data.at(i));
		}
		if(data.at(i).ch == 2) {
			ch2PESum += 1;
			pmtBHits.push_back(data.at(i));
		}
		
		// once we load our dataset, search forwards to find coincidences
		for(cur = i+1; cur < data.size(); cur++) { 
			
			// check the times of our two paired events. If the times are
			// not about the same, break. We don't have a coincidence! 
			if(data.at(cur).realtime - data.at(i).realtime > coincWindow*NANOSECOND) {
				break; 
			}
			// we only want coincidences in the dagger
			if(data.at(cur).ch != 1 && data.at(cur).ch != 2) { continue; }
			// count our coincidences in the same channel i 
			if(data.at(cur).ch == data.at(i).ch) {
				if(data.at(cur).ch == 1) {
					ch1PESum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == 2) {
					ch2PESum += 1;
					pmtBHits.push_back(data.at(cur));
				}
			}
			// count our coincidences in different channels
			if(data.at(cur).ch != data.at(i).ch) {				
				if(data.at(cur).ch == 1) {
					ch1PESum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == 2) {
					ch2PESum += 1;
					pmtBHits.push_back(data.at(cur));
				}
				// integrate the tail end. Add data to the sums of the 
				// two channels. 
				for(tailIt = cur+1; tailIt < data.size(); tailIt++) {
					if(data.at(tailIt).realtime - data.at(i).realtime > peSumWindow*NANOSECOND) {
						break;
					}
					if(data.at(tailIt).ch == 1) {
						ch1PESum += 1;
						pmtAHits.push_back(data.at(tailIt));
					}
					if(data.at(tailIt).ch == 2) {
						ch2PESum += 1;
						pmtBHits.push_back(data.at(tailIt));
					}
				}
				// Check to see if we found a neutron
				if(ch1PESum + ch2PESum > peSum) {
					coinc.push_back(data.at(i));
					coincIndices.push_back(i);
					pmtACoincHits.push_back(pmtAHits);
					pmtBCoincHits.push_back(pmtBHits);
					phsA.Fill(ch1PESum);
					phsB.Fill(ch2PESum);
					// deadtime correction 
					i = tailIt-1; 
					cur = tailIt-1;
					break;
				}
				break;
			}
		}
	}
	
	int numPhotonA = 0;
	int numPhotonB = 0;
	for(i = 1; i < (coinc.size()-1); i++) {
		if(coinc.at(i).realtime > 260 && coinc.at(i).realtime < 280) {
			if(((coinc.at(i+1).realtime-coinc.at(i).realtime) > 40000*NANOSECOND) && ((coinc.at(i).realtime-coinc.at(i-1).realtime) > 40000*NANOSECOND)) {
				numPhotonA = 0;
				numPhotonB = 0;
				cur = coincIndices.at(i);
				for(waveIt = coincIndices.at(i); waveIt < data.size(); waveIt++) {
					if(data.at(waveIt).realtime - data.at(cur).realtime > 40000*NANOSECOND) {
						printf("Data - %d,%d\n", numPhotonA, numPhotonB);
						break;
					}
					if(data.at(waveIt).ch == 1) {
						pmt1SummedWaveform.Fill((data.at(waveIt).realtime - data.at(cur).realtime)/NANOSECOND);
						allpmtAHits.push_back((data.at(waveIt).realtime - data.at(cur).realtime)/NANOSECOND);
						numPhotonA++;
					}
					if(data.at(waveIt).ch == 2) {
						pmt2SummedWaveform.Fill((data.at(waveIt).realtime - data.at(cur).realtime)/NANOSECOND);
						allpmtBHits.push_back((data.at(waveIt).realtime - data.at(cur).realtime)/NANOSECOND);
						numPhotonB++;
					}
				}
			}
		}
	}
	
	
}

/* Another coincidence timer. The difference between this one and the 
 * other type is that this one keeps track of the previous event. */
/*void Run::findMovingHists() {
	
	// check the data table and load it into the ROOT tree
	if(data.empty()) {
		this->readDataRoot();
	}
	if(data.empty()) {
		return;
	}

	// initialize iterators and variables
	int i;
	int cur = 0;
	int tailIt = 0;
	int waveIt = 0;
	
	// initialize variables for our PMT hits
	int ch1PESum;
	int ch2PESum;
	std::vector<input_t> pmtAHits;
	std::vector<input_t> pmtBHits;
	std::vector<long> coincIndices;

	// look at each entry of the data vector
	for(i = 0; i < data.size(); i++) {
		
		// clear previous buffers 
		ch1PESum = 0;
		ch2PESum = 0;
		pmtAHits.clear();
		pmtBHits.clear();
		
		input_t prevEvt;
		
		// for coincidence measurements we only want dagger hits
		if(data.at(i).ch != 1 && data.at(i).ch != 2) { continue; }
		
		// load info into our channel hit vectors
		if(data.at(i).ch == 1) {
			ch1PESum += 1;
			pmtAHits.push_back(data.at(i));
		}
		if(data.at(i).ch == 2) {
			ch2PESum += 1;
			pmtBHits.push_back(data.at(i));
		}
		
		// once we've loaded our dataset, search forwards to find coincidence
		for(cur = i+1; cur < data.size(); cur++) {
			
			// check the times of our two events. If the two are too far
			// apart, break because it's not a coincidence!
			if(data.at(cur).realtime - data.at(i).realtime > coincWindow*NANOSECOND) {
				break;
			}
			// only count coincidences in the dagger
			if(data.at(cur).ch != 1 && data.at(cur).ch != 2) { continue; }
			// count coincidences if we find a second event in channel i
			if(data.at(cur).ch == data.at(i).ch) {
				if(data.at(i).ch == 1) {
					ch1PESum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(i).ch == 2) {
					ch2PESum += 1;
					pmtBHits.push_back(data.at(cur));
				}
			}
			// count our coincidences on different channels
			if(data.at(cur).ch != data.at(i).ch) {
				if(data.at(cur).ch == 1) {
					ch1PESum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == 2) {
					ch2PESum += 1;
					pmtBHits.push_back(data.at(cur));
				}
				
				// save the previous event as a new data point
				prevEvt = data.at(cur);
				
				// integrate the tail end. Add counts into the right 
				// channel 
				for(tailIt = cur+1; tailIt < data.size(); tailIt++) {
					if((data.at(tailIt).realtime - prevEvt.realtime) > peSumWindow*NANOSECOND) {
						break;
					}
					if(data.at(tailIt).ch != 1 && data.at(tailIt).ch != 2) { continue; }
					if(data.at(tailIt).ch == 1) {
						ch1PESum += 1;
						pmtAHits.push_back(data.at(tailIt));
					}
					if(data.at(tailIt).ch == 2) {
						ch2PESum += 1;
						pmtBHits.push_back(data.at(tailIt));
					}
					prevEvt = data.at(tailIt);
				}
				// check to see if we've found a neutron!
				if(ch1PESum + ch2PESum >= peSum) {
					// fill our vectors
					coinc.push_back(data.at(i));
					coincIndices.push_back(i);
					pmtACoincHits.push_back(pmtAHits);
					pmtBCoincHits.push_back(pmtBHits);
					//phsA.Fill(ch1PESum);
					//phsB.Fill(ch2PESum);
					// put deadtime on counted neutrons
					i = tailIt-1; 
					cur = tailIt-1;
					break;
				}
				break;
			}
		}
	}
}
/* Removing code to make it easier to read
 * //				if(data.at(cur).ch == 1) {pmtAHits.push_back(data.at(cur));}
//				if(data.at(cur).ch == 2) {pmtBHits.push_back(data.at(cur));}
//				if(data.at(i).ch == 1) {pmtAHits.push_back(data.at(i));}
//				if(data.at(i).ch == 2) {pmtBHits.push_back(data.at(i));}
//				ch1PESum += 1;
//				ch2PESum += 1;
//if(ch1PESum >= peSum && ch2PESum >= peSum) { 
* //printf("Found %d coincidences\n", (int)coinc.size());

	//coincTiming->Write();
	//coincSpectrum->Write();
* //TH1F* coincTiming = new TH1F("coincTiming", "coincTiming", highWindow-lowWindow, lowWindow, highWindow); //Holds a timing spectrum of beta time - coincidence time
	//TH1F* coincSpectrum = new TH1F("coincSpectrum", "coincSpectrum", highEn-lowEn, lowEn, highEn); //Holds a PHS spectrum of gamma energies
	//Holds the calibration of energy (also masks out NaI for comparison)
	//printf("Finding Coincidences - Moving Winow\n");
//				if(data.at(cur).ch == 1) {pmtAHits.push_back(data.at(cur));}
//				if(data.at(cur).ch == 2) {pmtBHits.push_back(data.at(cur));}
//				if(data.at(i).ch == 1) {pmtAHits.push_back(data.at(i));}
//				if(data.at(i).ch == 2) {pmtBHits.push_back(data.at(i));}
//				ch1PESum += 1;
//				ch2PESum += 1;
	//if(ch1PESum >= peSum && ch2PESum >= peSum) { 
//printf("Found %d coincidences\n", (int)coinc.size());
	


	//coincTiming->Write();
	//coincSpectrum->Write();
*/
