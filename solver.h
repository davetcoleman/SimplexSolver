#include "utility.h"
#include <limits>

// Benchmarking
void runTests();

// Check if objective function is correct
void checkObjective(dictionary& d1, double answer);
	
// Run an occurance of linear solver
void solveLP(dictionary& d1);

// Run simplex
void simplex(dictionary& d1);

// Setup remainder of dictionary after it is loaded from file
void setup(dictionary& d1);

// Get the row index of the corresponding leaving vaiable
void getLeavingVar(dictionary& d1, int entering_var_index,
				   int& leaving_var_index, int& leaving_var_bound,
				   bool useBland, bool& doFlip);
	
// Get the row index of the next entering variable
void getEnteringVar(dictionary& d1, bool useBland, int& entering_var_index);
	 
// Check if problem is optimal
bool isOptimal(dictionary& d1);

// Check if dictionaty is feasible
bool isFeasible(dictionary& d1);

// Calculates the slack variables, ie the value of each row
void calculateSlack(dictionary& d1, bool checkObjValue);

// Do the actual pivot
void pivot(dictionary& d1, bool useBland);

// Dualize if initialization needed
void initialize(dictionary& d1);

// Combine original objective function to auxilary dictionary
void combineObjFunc(dictionary& d1, dictionary& s2);
