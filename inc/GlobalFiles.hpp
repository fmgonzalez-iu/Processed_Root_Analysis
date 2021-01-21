#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "stdio.h"

/*----------------------------------------------------------------------
 * Author: F. Gonzalez
 * 
 * This class is designed to hold paths to whatever pointers we might 
 * want to use for various running conditions. In particular, this is 
 * designed for like filling times or background rates. 
 * 
 * These are intended to be loaded initially through global environment
 * variables, though I'll leave it as an exercise to the reader to try 
 * and incorporate a functionality loading thing here.
 * 
 * I'm also going to (for now) exclude RUNDIR as one of these, since that
 * could be changed later.

 *--------------------------------------------------------------------*/


class GlobalFiles {
	// A class for loading files from global variables for e.g. traces
	// or filling time constants
		
	private:
	
		const char* fillLoc;  // MAD fill fit time const.
		const char* traceLoc; // Exp fill time const.
		const char* saveLoc;  // Save location for ROOT drawings
			
	public:
	
		GlobalFiles();
		~GlobalFiles();
	
};
