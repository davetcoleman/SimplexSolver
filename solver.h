#include "utility.h"
#include <limits>

// Get the row index of the corresponding leaving vaiable
int getLeavingVar(simplex s1, int entering_var_index);
	
// Get the row index of the next entering variable
int getEnteringVar(simplex s1);
	 
// Check if problem is optimal
bool isOptimal(simplex s1);

// Check if dictionaty is feasible
bool isFeasible(simplex s1);

// Return the bound that the nonbasic variable is currently resting on, or "active'
double getNonbasicVal(simplex s1, int index);
	
