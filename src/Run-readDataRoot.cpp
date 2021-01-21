#include "../inc/Run.hpp"
#include <algorithm>
#include "TLeaf.h"
#include "TBranch.h"
#include "TList.h"
#include "TKey.h"

#define NANOSECOND .000000001
#define CLKTONS 0.0000000008

/*	------------------------------------------------------------------------------------------------
	Author: Nathan B. Callahan (?)
	Editor: Frank M. Gonzalez
 
	This function takes the data vector and the file object and reads the data. First it discards
	header information, and then just slurps up the data into a struct which is appended to the end
	of the vector.
	------------------------------------------------------------------------------------------------	*/
	
int Run::numBits(uint32_t i)
{
     // Didn't come up with this, I assume it's magic. 
     i = i - ((i >> 1) & 0x55555555);
     i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
     return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

/* Load the data from ROOT into a format that C++ can use */
void Run::readDataRoot(const char* namecycle) {
	
	if (dataFile->IsZombie()) { // Make sure we didn't read a corrupted file
		return;
	}
	
	// initialize variables and ROOT tree 
	int numEntries;
	input_t event;
	event.deadtime = 16*NANOSECOND; // Hard-code in deadtime of 16 ns 
	// Technically, the deadtime for non-PMT events is variable. E.g. AC is 30-something ns.
	// This shouldn't matter, because rate dependence is most important for coincidences,
	// and monitor counts aren't really affected at our precision.
	event.rate = 0; // Also hard-code in rate of 0. This is used for random other stuff
	int i;
	TTree* rawData = NULL;
	
	// make sure we've initialized this properly
	if(this->exists() == false) {
		return;
	}
	
	// call the TList from the data file, to check the elements of our root tree 
	TList* list = dataFile->GetListOfKeys();
	if(list == NULL) {
		return;
	}	
	if(list->Contains(namecycle)) {	
		dataFile->GetObject(namecycle, rawData);
		// check if we got a tree from file
		if(rawData != NULL) {
			// load leaves from raw data tree
			TLeaf* timelf = rawData->FindLeaf("time");
			TLeaf* realtimelf = rawData->FindLeaf("realtime");
			TLeaf* chlf = rawData->FindLeaf("channel");
			TLeaf* taglf = rawData->FindLeaf("tag");
			
			// transfer the leaves to the C++ event
			timelf->SetAddress(&event.time);
			realtimelf->SetAddress(&event.realtime);
			chlf->SetAddress(&event.ch);
			taglf->SetAddress(&event.tag);
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
		}
		// sort the data, assuming we have data
		if(!data.empty()) {
			std::sort(data.begin(), data.end(), [](input_t x, input_t y)->bool{return (x.realtime < y.realtime);});
		}
	}
	return;
}
	// sort the data to a useful form 
	//if(!data.empty()) {
	//	std::sort(data.begin(), data.end(), [](input_t x, input_t y)->bool{return (x.realtime < y.realtime);});
	//}

	// close open root files to save memory
	//if(rawData) { delete rawData; }
	//dataFile->Close();
//	return;
//}

/* another function to read the ROOT data */
/*void Run::readDataRoot() {
	
	// initialize variables. dt is deadtime
	int numEntries;
	input_t event;
	int i;
	int dt;
	
	// initialize ROOT tree
	TTree* rawData = NULL;
	if(this->exists() == false) {
		return;
	}
	
	// load info from ROOT list
	TList* list = dataFile->GetListOfKeys();
	if(list == NULL) {return;}
	
	// preallocate namecycle memory
	char namecycle[32];
	// we can have multiple trees. We use the tree "tmcs_0" here 
	if(list->Contains("tmcs_1")) {
		sprintf(namecycle, "tmcs_1");	
		fprintf(stderr, "Reading tmcs_1 as main tree!!!\n");
		dataFile->GetObject(namecycle, rawData);
		
		// check if we got a tree from the file 
		if(rawData != NULL) {
			// import the leaves from our file
			TLeaf* timelf = rawData->FindLeaf("time");
			TLeaf* realtimelf = rawData->FindLeaf("realtime");
			TLeaf* chlf = rawData->FindLeaf("channel");
			TLeaf* taglf = rawData->FindLeaf("tag");
			// export the leaf data to a C++ event
			timelf->SetAddress(&event.time);
			realtimelf->SetAddress(&event.realtime);
			chlf->SetAddress(&event.ch);
			taglf->SetAddress(&event.tag);
		}
		// if no raw data return with empty vectors
		else {
			return;
		}

		// loop through the total entries and find their realtimes
		numEntries = rawData->GetEntries();
		for(i = 0; i < numEntries; i++) {
			rawData->GetEntry(i);
			event.realtime = ((double)event.time) * CLKTONS;
			data.push_back(event);
		}
	}

	// our second choice ROOT tree is mcs_events 
	else if(list->Contains("mcs_events")) {
		// we now have to load in the cycle # for the mcs_events tree
		TKey* key = (TKey*)list->First();
		if(key == NULL) {return;}
		int cycle = key->GetCycle();
		sprintf(namecycle, "mcs_events;%d", cycle);
		printf("Namecylce: %s\n", namecycle);

		// with our newly defined cycle name, we can now load the raw data
		dataFile->GetObject(namecycle, rawData);
	
		// check if we got a tree from the file
		if(rawData != NULL) {
			// load leaves and change them to events
			TLeaf* timelf = rawData->FindLeaf("time");
			TLeaf* realtimelf = rawData->FindLeaf("realtime");
			TLeaf* chlf = rawData->FindLeaf("ch");
			TLeaf* taglf = rawData->FindLeaf("tag");
			timelf->SetAddress(&event.time);
			realtimelf->SetAddress(&event.realtime);
			chlf->SetAddress(&event.ch);
			taglf->SetAddress(&event.tag);
		}
		// otherwise just return with empty vectors
		else {
			return;
		}
		
		// find the total number of entries
		numEntries = rawData->GetEntries();

		// loop through all the events 
		int flag = 0;
		for(i = 0; i < numEntries; i++) {
			rawData->GetEntry(i);
			// need to software-correct for multiple pulsing 
			if(event.ch == 5 && i > 0) {
				// find previous event on Ch. 5 
				for(dt = data.size()-1; dt >= 0; dt--) {
					if(data.at(dt).ch == 5) {
						break;
					}
				}
				// if the time between the most recent Ch. 5 evt is < DEADTIME, 
				// continue without putting in data . If dt is 0, then we 
				// had the first event. DEADTIME = 10 us
				if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 10000*NANOSECOND) {
						continue;
				}
			}
			// check if the event is in  channel 3. If so we need to find
			// how many tags we have on it. 
			if(event.ch == 3) {
				uint32_t tag = event.tag & (0x7800);
				int numTags = numBits(tag);
				// check for 3+ tag events
				if(numTags > 2) {
					printf("Found 3 tag event!\n");
					// find the previous Ch. 5 event
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 5 || data.at(dt).ch == 3) {
							break;
						}
					}
					// if we find a close event, it's probably one of the 
					// defining tag bit events.
					if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// map out previous bit
						tag = tag ^ data.at(dt).tag;
					}
					else {
						continue;
					}
				}
				// check for 0 tag events 
				else if(numTags == 0) {
					// find previous Ch. 5 event 
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 5 || data.at(dt).ch == 3) {
							break;
						}
					}
					// if we find a close event, it's probably the event that 
					// defines half of the tag bit
					if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// makes current tag the same as the previous event 
						tag = (1 << (data.at(dt).ch+5));
					}
				}
				// check for 2 tag events 
				else if(numTags == 2) {
					// find previous Ch. 3 event 
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 5 || data.at(dt).ch == 3) {
							break;
						}
					}
					// if this was a double followed by a double, then 
					// break them out and assign one channel to each. 
					if(dt >= 0 && numBits(data.at(dt).tag & (0x7800)) == 2 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						continue;
					}
					// if we find an event in close proximity, it's probably
					// the event which defines half of the tag bit
					else if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// map out previous bit 
						tag = tag ^ (1 << (data.at(dt).ch+5));
					}
					// if nothing else works, we can just loop around channels
					else {
						int t = 11;
						while(!(tag & (1<<t))) {
							t++;
						}
						event.ch = t - 5;
						data.push_back(event);
						tag = (tag ^ (1<<t));
					}
				}
				// use the tag to find which channel we're in
				switch(tag) {
					case (1 << 11) :
						event.ch = 6;
						break;
					case (1 << 12) :
						event.ch = 7;
						break;
					case (1 << 13) :
						event.ch = 8;
						break;
					case (1 << 14) :
						event.ch = 9;
						break;
					case 0:
						printf("Found 0 Tag event with no previous event within Gate Window!\n");
					default :
						break;
				}
			}
			// check what happens in channel 4 
			if(event.ch == 4) {
				uint32_t tag = event.tag & (0x600);
				

				int numTags = numBits(tag);
				if(numTags > 2) {
					printf("Found 3 tag event!\n");
					// check for the previous Ch. 5 event
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 9 || data.at(dt).ch == 4) {
							break;
						}
					}
					// If we find an event in close proximity, it's probably 
					// the event which defines half of the tag bit.
					if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// map out the previous bit 
						tag = tag ^ data.at(dt).tag;
					}
					else {
						continue;
					}
				}
				else if(numTags == 0) {
					// find the previous Ch. 5 event 
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 9 || data.at(dt).ch == 4) {
							break;
						}
					}
					// If we find an event in close proximity, it's probably 
					// the event which defines half of the tag bit.
					if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// make current tag same as previous event
						tag = (1 << (data.at(dt).ch-1));
					}
				}
				else if(numTags == 2) {
					// find previous Ch.3 event 
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 9 || data.at(dt).ch == 4) {
							break;
						}
					}
					// If this was a double followed by a double, then 
					// break them out and assign one channel to each 
					if(dt >= 0 && numBits(data.at(dt).tag & (0x600)) == 2 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						continue;
					}
					// If we find an event in close proximity, it's probably 
					// the event which defines half of the tag bit.
					else if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// map out previous bit 
						tag = tag ^ (1 << (data.at(dt).ch-1));
					}
					else {
						// if the event is outside, then scan and put into 
						// tag bits 
						int t = 9;
						while(!(tag & (1<<t))) {
							t++;
						}
						event.ch = t + 1;
						data.push_back(event);
						tag = (tag ^ (1<<t));
					}
				}
				// hardcode in channels for other tags
				switch(tag) {
					case (1 << 9) :
						event.ch = 10;
						break;
					case (1 << 10) :
						event.ch = 11;
						break;
					case 0:
						printf("Found 0 Tag event with no previous event within Gate Window!\n");
					default :
						break;
				}
			}
			data.push_back(event);
		}
		// initialize data iterators and initial data
		auto it = data.begin();
		auto backIt = data.begin();
		auto end = data.end();
		auto beg = data.begin();

		// go through the vector and impose a deadtime where appropriate,
		// depending on the channel. 
		for(it = data.begin(); it < end; it++) {
			// need software corrections for multiple pulsing 
			if((*it).ch == 9) {
				// find previous Ch. 5 event 
				for(backIt = it-1; backIt >= beg; backIt--) {
					if((*backIt).ch == 9) {
						break;
					}
				}
				// If the time between the most recent Ch.5 evt is < DEADTIME,
				// continue without putting in data. If dt is 0, then we
				// had the 1st evt.
				if(backIt >= beg && ((*it).realtime - (*backIt).realtime) < 10000*NANOSECOND) {
						(*backIt).ch=19;
						continue;
				}
			}
		}
	}
	
	// sort the data to a useful form 
	if(!data.empty()) {
		std::sort(data.begin(), data.end(), [](input_t x, input_t y)->bool{return (x.realtime < y.realtime);});
	}

	// close open root files to save memory
	if(rawData) { delete rawData; }
	dataFile->Close();
	return;
}*/

/* another function to read the ROOT data */
/*void Run::readDataEMS() {
	
	// initialize variables. dt is deadtime
	int numEntries;
	input_t event;
	int i;
	int dt;
	
	// initialize ROOT tree
	TTree* rawData = NULL;
	if(this->exists() == false) {
		return;
	}
	
	// load info from ROOT list
	TList* list = dataFile->GetListOfKeys();
	if(list == NULL) {return;}
	
	// preallocate namecycle memory
	char namecycle[32];
	// we can have multiple trees. We use the tree "tmcs_1" here 
	if(list->Contains("tmcs_1")) {
		sprintf(namecycle, "tmcs_1");	
		fprintf(stderr, "Reading tmcs_1 as main tree!!!\n");
		dataFile->GetObject(namecycle, rawData);
		
		// check if we got a tree from the file 
		if(rawData != NULL) {
			// import the leaves from our file
			TLeaf* timelf = rawData->FindLeaf("time");
			TLeaf* realtimelf = rawData->FindLeaf("realtime");
			TLeaf* chlf = rawData->FindLeaf("channel");
			TLeaf* taglf = rawData->FindLeaf("tag");
			// export the leaf data to a C++ event
			timelf->SetAddress(&event.time);
			realtimelf->SetAddress(&event.realtime);
			chlf->SetAddress(&event.ch);
			taglf->SetAddress(&event.tag);
		}
		// if no raw data return with empty vectors
		else {
			return;
		}

		// loop through the total entries and find their realtimes
		numEntries = rawData->GetEntries();
		for(i = 0; i < numEntries; i++) {
			rawData->GetEntry(i);
			event.realtime = ((double)event.time) * CLKTONS;
			data.push_back(event);
		}
	}

	// our second choice ROOT tree is mcs_events 
	else if(list->Contains("mcs_events")) {
		// we now have to load in the cycle # for the mcs_events tree
		TKey* key = (TKey*)list->First();
		if(key == NULL) {return;}
		int cycle = key->GetCycle();
		sprintf(namecycle, "mcs_events;%d", cycle);
		printf("Namecylce: %s\n", namecycle);

		// with our newly defined cycle name, we can now load the raw data
		dataFile->GetObject(namecycle, rawData);
	
		// check if we got a tree from the file
		if(rawData != NULL) {
			// load leaves and change them to events
			TLeaf* timelf = rawData->FindLeaf("time");
			TLeaf* realtimelf = rawData->FindLeaf("realtime");
			TLeaf* chlf = rawData->FindLeaf("ch");
			TLeaf* taglf = rawData->FindLeaf("tag");
			timelf->SetAddress(&event.time);
			realtimelf->SetAddress(&event.realtime);
			chlf->SetAddress(&event.ch);
			taglf->SetAddress(&event.tag);
		}
		// otherwise just return with empty vectors
		else {
			return;
		}
		
		// find the total number of entries
		numEntries = rawData->GetEntries();

		// loop through all the events 
		int flag = 0;
		for(i = 0; i < numEntries; i++) {
			rawData->GetEntry(i);
			// need to software-correct for multiple pulsing 
			if(event.ch == 5 && i > 0) {
				// find previous event on Ch. 5 
				for(dt = data.size()-1; dt >= 0; dt--) {
					if(data.at(dt).ch == 5) {
						break;
					}
				}
				// if the time between the most recent Ch. 5 evt is < DEADTIME, 
				// continue without putting in data . If dt is 0, then we 
				// had the first event. DEADTIME = 10 us
				if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 10000*NANOSECOND) {
						continue;
				}
			}
			// check if the event is in  channel 3. If so we need to find
			// how many tags we have on it. 
			if(event.ch == 3) {
				uint32_t tag = event.tag & (0x7800);
				int numTags = numBits(tag);
				// check for 3+ tag events
				if(numTags > 2) {
					printf("Found 3 tag event!\n");
					// find the previous Ch. 5 event
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 5 || data.at(dt).ch == 3) {
							break;
						}
					}
					// if we find a close event, it's probably one of the 
					// defining tag bit events.
					if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// map out previous bit
						tag = tag ^ data.at(dt).tag;
					}
					else {
						continue;
					}
				}
				// check for 0 tag events 
				else if(numTags == 0) {
					// find previous Ch. 5 event 
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 5 || data.at(dt).ch == 3) {
							break;
						}
					}
					// if we find a close event, it's probably the event that 
					// defines half of the tag bit
					if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// makes current tag the same as the previous event 
						tag = (1 << (data.at(dt).ch+5));
					}
				}
				// check for 2 tag events 
				else if(numTags == 2) {
					// find previous Ch. 3 event 
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 5 || data.at(dt).ch == 3) {
							break;
						}
					}
					// if this was a double followed by a double, then 
					// break them out and assign one channel to each. 
					if(dt >= 0 && numBits(data.at(dt).tag & (0x7800)) == 2 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						continue;
					}
					// if we find an event in close proximity, it's probably
					// the event which defines half of the tag bit
					else if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// map out previous bit 
						tag = tag ^ (1 << (data.at(dt).ch+5));
					}
					// if nothing else works, we can just loop around channels
					else {
						int t = 11;
						while(!(tag & (1<<t))) {
							t++;
						}
						event.ch = t - 5;
						data.push_back(event);
						tag = (tag ^ (1<<t));
					}
				}
				// use the tag to find which channel we're in
				switch(tag) {
					case (1 << 11) :
						event.ch = 6;
						break;
					case (1 << 12) :
						event.ch = 7;
						break;
					case (1 << 13) :
						event.ch = 8;
						break;
					case (1 << 14) :
						event.ch = 9;
						break;
					case 0:
						printf("Found 0 Tag event with no previous event within Gate Window!\n");
					default :
						break;
				}
			}
			// check what happens in channel 4 
			if(event.ch == 4) {
				uint32_t tag = event.tag & (0x600);
				

				int numTags = numBits(tag);
				if(numTags > 2) {
					printf("Found 3 tag event!\n");
					// check for the previous Ch. 5 event
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 9 || data.at(dt).ch == 4) {
							break;
						}
					}
					// If we find an event in close proximity, it's probably 
					// the event which defines half of the tag bit.
					if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// map out the previous bit 
						tag = tag ^ data.at(dt).tag;
					}
					else {
						continue;
					}
				}
				else if(numTags == 0) {
					// find the previous Ch. 5 event 
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 9 || data.at(dt).ch == 4) {
							break;
						}
					}
					// If we find an event in close proximity, it's probably 
					// the event which defines half of the tag bit.
					if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// make current tag same as previous event
						tag = (1 << (data.at(dt).ch-1));
					}
				}
				else if(numTags == 2) {
					// find previous Ch.3 event 
					for(dt = data.size()-1; dt >= 0; dt--) {
						if(data.at(dt).ch > 9 || data.at(dt).ch == 4) {
							break;
						}
					}
					// If this was a double followed by a double, then 
					// break them out and assign one channel to each 
					if(dt >= 0 && numBits(data.at(dt).tag & (0x600)) == 2 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						continue;
					}
					// If we find an event in close proximity, it's probably 
					// the event which defines half of the tag bit.
					else if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 1000*NANOSECOND) {
						// map out previous bit 
						tag = tag ^ (1 << (data.at(dt).ch-1));
					}
					else {
						// if the event is outside, then scan and put into 
						// tag bits 
						int t = 9;
						while(!(tag & (1<<t))) {
							t++;
						}
						event.ch = t + 1;
						data.push_back(event);
						tag = (tag ^ (1<<t));
					}
				}
				// hardcode in channels for other tags
				switch(tag) {
					case (1 << 9) :
						event.ch = 10;
						break;
					case (1 << 10) :
						event.ch = 11;
						break;
					case 0:
						printf("Found 0 Tag event with no previous event within Gate Window!\n");
					default :
						break;
				}
			}
			data.push_back(event);
		}
		// initialize data iterators and initial data
		auto it = data.begin();
		auto backIt = data.begin();
		auto end = data.end();
		auto beg = data.begin();

		// go through the vector and impose a deadtime where appropriate,
		// depending on the channel. 
		for(it = data.begin(); it < end; it++) {
			// need software corrections for multiple pulsing 
			if((*it).ch == 9) {
				// find previous Ch. 5 event 
				for(backIt = it-1; backIt >= beg; backIt--) {
					if((*backIt).ch == 9) {
						break;
					}
				}
				// If the time between the most recent Ch.5 evt is < DEADTIME,
				// continue without putting in data. If dt is 0, then we
				// had the 1st evt.
				if(backIt >= beg && ((*it).realtime - (*backIt).realtime) < 10000*NANOSECOND) {
						(*backIt).ch=19;
						continue;
				}
			}
		}
	}
	
	// sort the data to a useful form 
	if(!data.empty()) {
		std::sort(data.begin(), data.end(), [](input_t x, input_t y)->bool{return (x.realtime < y.realtime);});
	}

	// close open root files to save memory
	if(rawData) { delete rawData; }
	dataFile->Close();
	return;
}

/* Removed (commented) code for cleanliness):
 * //int numKeys = dataFile->GetNkeys();
			//TBranch* br = rawData->GetBranch("events");
//	event.realtime = ((double)event.time) * CLKTONS;
//int numKeys = dataFile->GetNkeys();
* //	if(list->Contains("tmcs_0")) {
//		sprintf(namecycle, "tmcs_0");	
//TBranch* br = rawData->GetBranch("events");

//dataFile->GetObject("mcs_events", rawData);
//TBranch* br = rawData->GetBranch("events");
//numMonInEntry = 0;
* /*if(event.realtime > 1000) {flag = 1;}
			if(flag==1){printf("%f\n", event.realtime);}*/
//if(dt >= 0 && (event.realtime - data.at(dt).realtime) < 10000*NANOSECOND) { //If the time between the most recent Ch.5 evt is < DEADTIME, continue without putting in data. If dt is 0, then we had the 1st evt.
//printf("%d | %d %d %d %d | %4.9f\n", numBits(tag), (tag & (1 << 11))>>11, (tag & (1 << 12))>>12, (tag & (1 << 13))>>13, (tag & (1 << 14))>>14, event.realtime);
//printf("Found 0 tag event!\n");
//printf("Found event for 0-Mux-Bit within %f\n", (event.realtime - data.at(dt).realtime)/NANOSECOND);
//printf("%d | %d %d %d %d | %4.9f\n", numBits(tag), (tag & (1 << 11))>>11, (tag & (1 << 12))>>12, (tag & (1 << 13))>>13, (tag & (1 << 14))>>14, event.realtime);
//printf("Found 0 tag event!\n");
//printf("Found event for 0-Mux-Bit within %f\n", (event.realtime - data.at(dt).realtime)/NANOSECOND);
