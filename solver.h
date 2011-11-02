#include "utility.h"
#include <limits>

// Unit Testing?
void runTests();

// Compare two matricies
bool dictionaryIsEqual(dictionary d1, dictionary d);

// Run an occurance of linear solver
dictionary solveLP(dictionary s1);

// Run simplex
dictionary simplex(dictionary s1);

// Setup remainder of dictionary after it is loaded from file
dictionary setup(dictionary s1);

// Get the row index of the corresponding leaving vaiable
void getLeavingVar(dictionary s1, int entering_var_index, int& leaving_var_index, int& leaving_var_bound);
	
// Get the row index of the next entering variable
int getEnteringVar(dictionary s1);
	 
// Check if problem is optimal
bool isOptimal(dictionary s1);

// Check if dictionaty is feasible
bool isFeasible(dictionary s1);

// Calculates the slack variables, ie the value of each row
dictionary calculateSlack(dictionary s1);

// Do the actual pivot
dictionary pivot(dictionary s1);

// Dualize if initialization needed
dictionary initialize(dictionary s1);

// Combine original objective function to auxilary dictionary
dictionary combineObjFunc(dictionary s1, dictionary s2);
