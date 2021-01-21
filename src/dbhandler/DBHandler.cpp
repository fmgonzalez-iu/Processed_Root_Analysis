#include "../inc/DBHandler.hpp"

/*	------------------------------------------------------------------------------------------------
	Author: Nathan B. Callahan
	Editor: Frank M. Gonzalez
	
	This file contains the constructor and destructor for the DBHandler class. It sets the members
	to the given values and calls the private method getRuns() which queries the MySQL database and
	prepares the list of runs for iteration.
	------------------------------------------------------------------------------------------------	*/

/* database constructor */
DBHandler::DBHandler(const char* sqlQuery, int coincWindow, int peSumWindow, int peSum, int coincMode) {
	query = strdup(sqlQuery);
	this->coincWindow = coincWindow;
	this->peSumWindow = peSumWindow;
	this->peSum = peSum;
	this->coincMode = coincMode;
	this->getRuns();
}

/* database destructor */
DBHandler::~DBHandler() {
	delete query;
}

/* get runs for iterations */
std::vector<double> DBHandler::getXs() {
	return xs;
}
