#include "../inc/GlobalFiles.hpp"

GlobalFiles::GlobalFiles() {
	
	this->fillLoc = std::getenv("FILL_LOC");
	this->traceLoc = std::getenv("TRACE_LOC");
	this->saveLoc = std::getenv("SAVE_LOC");
	
}

GlobalFiles::~GlobalFiles() {
	
	delete fillLoc;
	delete traceLoc;
	delete saveLoc;
	
}
