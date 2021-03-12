#include "../inc/Run.hpp"

/*----------------------------------------------------------------------
	Author: Nathan B. Callahan (?)
	Editor: Frank M. Gonzalez
	
These functions transform our run data into histograms. This is useful 
because it lets us play around with the data in ROOT, and import/export
the data files easier.	
-----------------------------------------------------------------------*/
//----------------------------------------------------------------------
/* Make vectors. These extract data from our .ROOT file 
 * This is useful for outputtting values and getting sums */
//----------------------------------------------------------------------
std::vector<counting_t> Run::getCountingAvgs(int nvar, const std::function <coinc_t (coinc_t)>& expr,
						const std::function <bool (coinc_t)>& selection) 
{
	std::vector<coinc_t> filtered;
	std::vector<coinc_t> transformed;
	
	if (coinc.empty()) { // Load counts if we haven't already.
		this->getCoincCounts(expr,selection);
	}
	if (nvar < 1) {
		printf("Warning! nVar for counting averages must be greater than 1!\n");
		nvar = 1;
	}
	
	// Filter the coincidence counts ala lambda algebra (like normal coincidences)
	std::copy_if(coinc.begin(), coinc.end(), std::back_inserter(filtered), selection);
	std::transform(filtered.begin(), filtered.end(), std::back_inserter(transformed), expr);

	// Set the initial average as just raw coincidence counts
	counting_t avg;
	avg.dt = std::accumulate(transformed.begin(), transformed.end(), 0.0, [](double t, coinc_t x)->double{return t + (x.length);}) / (double)transformed.size();// Average length of coincidence
	avg.N  = std::accumulate(transformed.begin(), transformed.end(), 0.0, [](double N, coinc_t x)->double{return N + (double)(x.pmt1 + x.pmt2);}) / (double)transformed.size(); // Average number of photons in coincidences

	std::vector<counting_t> averages;
	averages.push_back(avg); // This is of all coincidences -- average N and dt
	
	// Now we need to loop through nVar;
	// First, let's get timing and pe parameters.
	std::vector<unsigned long> tThresh;
	for (int ii = 1; ii < nvar+1; ii++) { 
		tThresh.push_back((unsigned long)((this->promptWindow * 0.8) * (2*ii))); // Integer multiples of promptWindow.
	}
	std::vector<int> peThresh;
	for (int ii = 1; ii < nvar+1; ii++) {
		peThresh.push_back((int)(this->peSum * (2*ii))); // Integer multiples of peSum
	}
	
	// Now we calculate the averages for each of these
	for (size_t ii = 0; ii < nvar; ii++) {	
		counting_t timing; // timing based average -- fixed time, non-fixed PE
		timing.dt = (double)(tThresh.at(ii) * NANOSECOND);
		timing.N  = 0.;
		counting_t photon; // PE based average -- fixed sum, non-fixed timing
		photon.dt = 0.;
		photon.N  = (double)peThresh.at(ii);
		printf("ii=%lu,dt=%e,N=%f\n",ii,(double)(tThresh.at(ii)*NANOSECOND),(double)peThresh.at(ii));
		double peCounts = 0.;
		double tCounts = 0.;
		for (auto cc = transformed.begin(); cc < transformed.end(); cc++) {
			// Get the counts from combined PMTs
			std::vector<input_t> pmt = this->pmtCCoincHits.at((*cc).index);
			/*std::vector<input_t> pmt = this->getCounts(
					[](input_t x)->input_t{return x;},
				[this,cc](input_t x)->bool{return ((x.ch == this->coinChan1 || x.ch == this->coinChan2)
													&& x.realtime >= (*cc).realtime && (x.realtime-(*cc).realtime) < 0.001);}
				);*/
			//double timing = 0.;
			unsigned long tIni = (*cc).time;
					
			//std::vector<input_t> pmt1 = this->pmtACoincHits.at((*cc).index);
			//std::vector<input_t> pmt2 = this->pmtBCoincHits.at((*cc).index);
			// Combine PMT 1 and PMT 2 into a single vector
			//std::vector<input_t> pmtC = pmt1;
			//pmtC.insert(pmtC, pmt2.begin(), pmt2.end());
			//pmtC.sort(
			
			// increment timing N
			timing.N = std::accumulate(pmt.begin(), pmt.end(), timing.N, 
							[tThresh, tIni, ii](double p, input_t x)->double{
								return p + ((x.time-tIni) < tThresh.at(ii) ? 1 : 0);});
			/*printf("accumulate=%d\n",std::accumulate(pmt.begin(), pmt.end(), 0, 
							[tThresh, ii](double p, input_t x)->double{
								return p + (x.time < tThresh.at(ii) ? 1 : 0);}));
			*/
			//timing.N = std::accumulate(pmt1.begin(), pmt1.end(), timing.N, 
							//[tThresh, ii](double p, input_t x)->double{
								//return p + (x.time < tThresh.at(ii) ? 1 : 0);});
			//timing.N = std::accumulate(pmt2.begin(), pmt2.end(), timing.N, 
							//[tThresh, ii](double p, input_t x)->double{
								//return p + (x.time < tThresh.at(ii) ? 1 : 0);});
			// and add PMT sum info
			if (pmt.size() > peThresh.at(ii)) {
				photon.dt += (pmt.at((size_t)peThresh.at(ii)).realtime - (*cc).realtime) + (this->peSumWindow*0.8*NANOSECOND);
			//	peCounts += 1;
			} else { 
				photon.dt += (pmt.back().realtime - (*cc).realtime) + (this->peSumWindow*0.8*NANOSECOND);
			}
			
			/*if (pmt1.size() + pmt2.size() > peThresh.at(ii)) {
				if (pmt1.at((size_t)peThresh.at(ii)).time > pmt2.at((size_t)peThresh.at(ii)).time) {
					photon.dt += pmt1.at((size_t)peThresh.at(ii)).realtime+(this->peSumWindow*0.8*NANOSECOND);
				} else {
					photon.dt += pmt2.at((size_t)peThresh.at(ii)).realtime+(this->peSumWindow*0.8*NANOSECOND);
				}
			} else {
				if (pmt1.back().time > pmt2.back().time) {
					photon.dt += pmt1.back().realtime+(this->peSumWindow*0.8*NANOSECOND);
				} else {
					photon.dt += pmt2.back().realtime+(this->peSumWindow*0.8*NANOSECOND);
				}
			}*/
		}
		timing.N  /= (double)(transformed.size()); // And take averages
		photon.dt /= (double)(transformed.size());
		averages.push_back(timing);
		averages.push_back(photon);
		// DEBUG
		printf("T: %f, %e\n",timing.N, timing.dt);
		printf("PE: %f, %e\n",photon.N, photon.dt);
		
	}
	return averages;
	
}						

/* Find the total number of coincidence counts in our data set (using coinc inputs) */
std::vector<coinc_t> Run::getCoincCounts(const std::function <coinc_t (coinc_t)>& expr,
					 const std::function <bool (coinc_t)>& selection)
{
	// load our data vectors 
	std::vector<coinc_t> filtered;
	std::vector<coinc_t> transformed;

	// These are defined in Run-findCoincidence.cpp.
	// TODO: Simplify getCoincCounts switch with findCoinc (but fiddling around with internal params	
	if (coinc.empty()) { // Assuming we haven't already loaded run-coinc
		if ((coincMode == 1) || (coincMode == 2)) { // Fixed
			this->findCoincidenceFixed();
		} else if ((coincMode == 3) || (coincMode == 4)) { // Telescoping
			this->findCoincidenceMoving();
		} else if ((coincMode == 5) || (coincMode == 6)) { // Telescoping, don't care which PMT.
			this->findCoincidenceNoPMT();
		} else if ((coincMode == 7) || (coincMode == 8)) { // Telescoping, but with additional singles deadtime
			this->findCoincidenceFixedTele();
		} else if ((coincMode == 9) || (coincMode == 10)) { // Telescoping, filter out high-PE coincidences
			this->findCoincidenceFixedTele();
		} else {
			printf("ERROR! Coincidence mode not defined in Run-makeHists.cpp!");
			return transformed; // Return nothing
		}
	}
		//if (coincMode <= 4) { // Normal running mode
			//if (coincMode % 2 == 0) {
				//this->findCoincidenceMoving();
			//} else {
				//this->findCoincidenceFixed();
			//}
		//} else if ((coincMode == 5) || (coincMode == 6)) {
			//this->findCoincidenceNoPMT(); // For no PMT
		//}
	//}
	
	// copy and transform our data sets to find the total amount of counts in the ROI
	std::copy_if(coinc.begin(), coinc.end(), std::back_inserter(filtered), selection);
	std::transform(filtered.begin(), filtered.end(), std::back_inserter(transformed), expr);
	printf("%lu\n",filtered.size());
	// Find the instantaneous rate of coincidences. (Maybe more accurate deadtime?)
	size_t ii = 0; // Iterator (where to start rate counting)
	for (size_t jj = 0; jj < transformed.size(); jj++) {
		
		// Want to find how many events are within 0.1s of this event
		double dt = 0.1; // 0.1 second averaging
		while ((transformed.at(jj).realtime - transformed.at(ii).realtime) >=  dt) { 
			ii+=1; 
		}
		if (ii > jj) { ii = jj; } // Can't have a negative size
		double cts = (double)(jj - ii); // How many counts in the window?
		
		transformed.at(jj).rate = cts/dt; // And set the rate to cts/dt.	
	}
	return transformed;
}

// Get AC/dagger coincidence counts
std::vector<coinc_t> Run::getCoincCountsAC(const std::function <coinc_t (coinc_t)>& expr,
				  const std::function <bool (coinc_t)>& selection)
{
	// load our data vectors 
	std::vector<coinc_t> filtered;
	std::vector<coinc_t> transformed;

	// find coincidences from our coincidence mode -- sketchily coded
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findCoincidenceAC();
		} 
		else if(coincMode == 2) {
			this->findAntiCoincidenceAC();
		}
	}
	
	// copy and transform our data sets to find the total amount of counts
	std::copy_if(coinc.begin(), coinc.end(), std::back_inserter(filtered), selection);
	std::transform(filtered.begin(), filtered.end(), std::back_inserter(transformed), expr);
	
	return transformed;
}

/* Find the individual counts in our data set (sorted, set by detector) */
std::vector<input_t> Run::getCounts(const std::function <input_t (input_t)>& expr,
				  const std::function <bool (input_t)>& selection)
{
	// load our data vectors
	std::vector<input_t> filtered;
	std::vector<input_t> transformed;
	
	// load our root tree 
	if(data.empty()) {
		//this->readDataRoot();
		input_t blank;
		transformed.push_back(blank);
		return transformed;
	}
		
	// copy and transform our data sets to find the total amount of counts
	std::copy_if(data.begin(), data.end(), std::back_inserter(filtered), selection);
	std::transform(filtered.begin(), filtered.end(), std::back_inserter(transformed), expr);
	
	// Calculate the instantaneous rate for each photon
	size_t ii = 0; // Iterator (where to start rate counting)
	for (size_t jj = 0; jj < transformed.size(); jj++) {
		
		// Want to find how many events are within 0.1s of this event
		double dt = 0.1; // 0.1 second averaging
		while ((transformed.at(jj).realtime - transformed.at(ii).realtime) >=  dt) { 
			ii+=1; 
		}
		if (ii > jj) { ii = jj; } // Can't have a negative size
		double cts = (double)(jj - ii); // How many counts in the window?
		
		transformed.at(jj).rate = cts/dt; // And set the rate to cts/dt.	
	}
		
	return transformed;
}

/* This function produces vectors of our photon's traces. */
/*std::vector<double> Run::getPhotonTracesVect(int pmt,
											 const std::function <double (std::vector<input_t>)>& expr,
											 const std::function <bool (std::vector<input_t>)>& selection)
{
	
	// initialize the photon vector, which requires some data from the 
	// coincidence hits vector. Also initialize the data vectors.
	std::vector<std::vector<input_t> > *photonVect = pmt == 1 ? &pmtACoincHits : &pmtBCoincHits;
	std::vector<std::vector<input_t> > filtered;
	std::vector<double> mapped;

	// load data from ROOT depending on what the coincidence mode is
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
	
	// see if the photon vector exists. If it doesn't, end. 
	if(photonVect->empty()) {
		return mapped;
	}
	
	// copy the photon vector to pick a particular selection from this
	std::copy_if(photonVect->begin(), photonVect->end(), std::back_inserter(filtered), selection);
	
	// close if we picked the wrong dataset 
	if(filtered.empty()) {
		return mapped;
	}
	
	// load a transform vector to modify our dataset and load the photon
	// traces
	std::transform(filtered.begin(), filtered.end(), std::back_inserter(mapped), expr);
	
	return mapped;
}*/

//----------------------------------------------------------------------
/* Make Root Histograms. This is useful for plotting and visualizations */
//----------------------------------------------------------------------
/* This function takes two inputs the functional expression and the 
 * selected data, and produces a histogram. */
/*TH1D Run::getCoincHist(const std::function <double (input_t)>& expr,
					   const std::function <bool (input_t)>& selection)
{
	
	// check to load our ROOT tree, depending on whether we are in singles
	// or doubles mode. 
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
	
	// close function if we select an empty set 
	if(coinc.empty()) {
		TH1D hist("Empty_histo", "Empty_histo", 10, 0, 10);
		return hist;
	}
	
	// initialize our data vectors (in a unique input_t class) 
	std::vector<input_t> filtered;
	std::copy_if(coinc.begin(), coinc.end(), std::back_inserter(filtered), selection);
	
	// close function if we select an empty set 
	if(filtered.empty()) {
		TH1D hist("Empty_histo", "Empty_histo", 10, 0, 10);
		return hist;
	}
	
	// load the min and max from our filtered dataset, with bounds set by 
	// some functional expressions on arbitrary x and y 
	input_t min = *std::min_element(filtered.begin(), filtered.end(), [expr](input_t x, input_t y)->bool{return(expr(x) < expr(y));});
	input_t max = *std::max_element(filtered.begin(), filtered.end(), [expr](input_t x, input_t y)->bool{return(expr(x) < expr(y));});
	
	// if we haven't closed out yet, create a histogram 
	TH1D hist("histo", "histo", ceil(expr(max))-floor(expr(min)), floor(expr(min)), ceil(expr(max)));
	
	// loop our histogram through our whole analysis code, and modify
	// histogram in our coincidence runs. 
	std::for_each(filtered.begin(), filtered.end(), [&hist, expr](input_t event)->void{hist.Fill(expr(event));});
	
	// output is our coincidence histogram 
	return hist;
}*/

/* This function takes the same inputs as our other coincidence hist
 * but doesn't actually require coincidences. Instead it just produces
 * an arbitrary hist. */
/*TH1D Run::getHist(const std::function <double (input_t)>& expr,
				  const std::function <bool (input_t)>& selection)
{
	// load our ROOT tree into the data 
	if(data.empty()) {
		this->readDataRoot();
	}
	
	// make sure we've actually picked a data set with data. If not, then
	// close with an empty histogram
	if(data.empty()) {
		TH1D hist("Empty_histo", "Empty_histo", 0, 0, 0);
		return hist;
	}
	
	// initialize data vectors 
	std::vector<input_t> filtered;
	std::copy_if(data.begin(), data.end(), std::back_inserter(filtered), selection);
	
	// check to make sure our data vectors actually exist. 
	if(filtered.empty()) {
		TH1D hist("Empty_histo", "Empty_histo", 0, 0, 0);
		return hist;
	}
	
	// load the min and max from our filtered data set, defined through 
	// some functional we're acting on with x and y 
	input_t min = *std::min_element(filtered.begin(), filtered.end(), [expr](input_t x, input_t y)->bool{return(expr(x) < expr(y));});
	input_t max = *std::max_element(filtered.begin(), filtered.end(), [expr](input_t x, input_t y)->bool{return(expr(x) < expr(y));});
	
	// initialize our histogram 
	TH1D hist("histo", "histo", ceil(expr(max))-floor(expr(min)), floor(expr(min)), ceil(expr(max)));
	
	// loop our histogram through the analyzer software. 
	std::for_each(filtered.begin(), filtered.end(), [&hist, expr](input_t event)->void{hist.Fill(expr(event));});
	
	// output our final histogram 
	return hist;
}*/

/* This histogram only takes the start and endtimes of a dataset, and 
 * uses these to produce a deadtime histogram. */
/*TH1D Run::getDeadtimeHist(double start, double end) {
	
	// initialize histogram with input ranges 
	TH1D deadTimeHist("deadTime", "deadTime", ceil(end)-floor(start), floor(start), ceil(end));
	
	// load in coincidence data from ROOT 
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
	
	// if the coincidence dataset is empty, return out hist with no hits 
	if(pmtACoincHits.empty()) {
		return deadTimeHist;
	}
	
	// initialize data variables. Assume the coincidence hits vectors 
	// were filled together and thus have the same number of entries 
	auto pmtAit = pmtACoincHits.begin();
	auto pmtBit = pmtBCoincHits.begin();
	std::vector<input_t> pmtAwaveform;
	std::vector<input_t> pmtBwaveform;
	double firstCountTime;
	double lastCountTime;
	double deadTime;
	
	// loop through the coincidence hits on pmtA. If A and B have the 
	// same length this will still work, we increment B at the end. 
	for(pmtAit = pmtACoincHits.begin(); pmtAit < pmtACoincHits.end(); pmtAit++) {
		
		// load our waveforms by calling the location of our iterators 
		pmtAwaveform = *pmtAit;
		pmtBwaveform = *pmtBit;

		// take the first and last count times for both A waveforms and B
		// waveforms. Take the minimum and maximum of both of these
		firstCountTime = std::min(pmtAwaveform.front().realtime, pmtBwaveform.front().realtime);
		lastCountTime = std::max(pmtAwaveform.back().realtime, pmtBwaveform.back().realtime);
		
		// check to make sure the count time is in the right bounds 
		if(firstCountTime > start && firstCountTime < end) {
			
			// increment the dead time histogram by adding in the difference
			// between the initial and final counting times
			deadTime = lastCountTime-firstCountTime;
			deadTimeHist.Fill(firstCountTime, deadTime);		
		}
		pmtBit++;
	}
	return deadTimeHist;
}*/
/* A function that loads histograms */
/*TH1D Run::getHistIterator(const std::function <double (std::vector<input_t>::iterator,
													   std::vector<input_t>::iterator,
													   std::vector<input_t>::iterator)>& expr,
						  const std::function <bool (std::vector<input_t>::iterator,
													   std::vector<input_t>::iterator,
													   std::vector<input_t>::iterator)>& selection)
{
	// check the data and read into ROOT 
	if(data.empty()) {
		this->readDataRoot();
	}
	
	// close if we accidentally load a bad dataset
	if(data.empty()) {
		TH1D hist("Empty_histo", "Empty_histo", 10, 0, 10);
		return hist;
	}
	
	// initialize data vectors and iterators 
	std::vector<double> points;
	std::vector<input_t>::iterator it;
	
	// loop through our data points and select the ones that we want in
	// our histogram 
	for(it = data.begin(); it < data.end(); it++) {
		if(selection(it, data.begin(), data.end()) == true) {
			points.push_back(expr(it, data.begin(), data.end()));
		}
	}
	
	// initialize min and max  
	double min;
	double max;
	
	// load minimum and maximum elements if we have data. If we don't 
	// have good data, then choose zero and one
	if(points.size() > 0) {
		min = *std::min_element(points.begin(), points.end());
		max = *std::max_element(points.begin(), points.end());
	}
	else {
		min = 0;
		max = 1;
	}

	// create our histogram 
	TH1D hist("histo", "histo", 1000, min, max);
	printf("min: %f; max: %f\n", min, max);
	
	// run our produced histograms through the analyzer
	if(points.size() > 0) {
		std::for_each(points.begin(), points.end(), [&hist](double point)->void{hist.Fill(point);});
	}
	return hist;
}*/

/*-----------------------------------------------------------------------------------------------
 * Extra code goes here
 * //printf("Count time, length: %e, %e\n", firstCountTime, lastCountTime-firstCountTime);
 * //deadTime = lastCountTime-firstCountTime > 500e-9 ? lastCountTime-firstCountTime : 500e-9;
 * //printf("Deadtime: %e\n", lastCountTime-firstCountTime);
 * /*std::vector<double> Run::getDeadtimeVect(double start, double end) {
	std::vector<double> dtVect;
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findcoincidenceMoving();
		}
	}
	
	if(pmtACoincHits.empty()) {
		return dtVect;
	}
	
	auto pmtAit = pmtACoincHits.begin();
	auto pmtBit = pmtBCoincHits.begin(); //We'll assume that these vectors were filled together and therefore have same # of entries
	std::vector<input_t> pmtAwaveform;
	std::vector<input_t> pmtBwaveform;
	double firstCountTime;
	double lastCountTime;
	double deadTime;
	
	for(pmtAit = pmtACoincHits.begin(); pmtAit < pmtACoincHits.end(); pmtAit++) {
		pmtAwaveform = *pmtAit;
		pmtBwaveform = *pmtBit;

		firstCountTime = std::min(pmtAwaveform.front().realtime, pmtBwaveform.front().realtime);
		lastCountTime = std::max(pmtAwaveform.back().realtime, pmtBwaveform.back().realtime);
		//printf("Count time, length: %e, %e\n", firstCountTime, lastCountTime-firstCountTime);
		
		if(firstCountTime > start && firstCountTime < end) {
			//deadTime = lastCountTime-firstCountTime > 500e-9 ? lastCountTime-firstCountTime : 500e-9;
			deadTime = lastCountTime-firstCountTime;
			deadTimeHist.Fill(firstCountTime, deadTime);
			//printf("Deadtime: %e\n", lastCountTime-firstCountTime);
		}
		pmtBit++;
	}
	
	return deadTimeHist;
}*/
			/*if(expr(it, data.begin(), data.end()) > 1.0) {
				printf("expr: %f", expr(it, data.begin(), data.end()));
			}*/
