#include "../inc/DBHandler.hpp"

/*	------------------------------------------------------------------------------------------------
	Author: Nathan B. Callahan
	Editor: Frank M. Gonzalez
	
	This file contains the method that queries the database and loads the runs.
	It works through a MySQL server, which allows for "easy" database management.
	------------------------------------------------------------------------------------------------	*/

void DBHandler::getRuns() {
	// connect to server
	TMySQLServer* serv = new TMySQLServer("mysql://localhost", "root", "iucf1234");
	TSQLResult* res = serv->Query(query);
	// check if we got 2 columns of data (runs and xs)
	if(res != NULL && res->GetFieldCount() != 2) {
		fprintf(stderr, "Error! Result of Mysql query did not have 2 columns! Stopping analysis\n\n*****The query MUST contain the run numbers and x values to plot (in order)*****\n\n");
		exit(1);
	}
	
	if(res == NULL) {
		fprintf(stderr, "Error! Got a NULL result. Stopping analysis\n");
		exit(1);
	}
	// gets how many rows were retrieved
	int num = res->GetRowCount();
	
	TSQLRow* row = NULL;
	
	int i;
	
	// loop through all rows 
	for(i = 0; i < num; i++) {
		row = res->Next();
		if(row == NULL) {
			fprintf(stderr, "Error! Got a NULL row when one should exist. Stopping analysis\n");
			exit(1);
		}
		
		// push onto runs and xs from the 0th and 1st columns of the row
		runs.push_back(atoi(row->GetField(0)));
		runBodies.push_back(row->GetField(1));
	}
	serv->Close();
}
