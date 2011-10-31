#include "utility.h"
#include <limits>

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

// Return the bound that the nonbasic variable is currently resting on, or "active'
double getNonbasicVal(dictionary s1, int index);

// Do the actual pivot
dictionary pivot(dictionary s1);

