// LP Solver for General Simplex Problems
// CSCI Linear Programming - Programming Assignment 1
// Dave Coleman | david.t.coleman@colorado.edu
// 10/20/2011
//
// Solve a simple linear program

//----------------------------------------------------------
// Directives:
//----------------------------------------------------------
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "armadillo"
#include <sstream>
#include <string>
#include "solver.h"
#include "assert.h"


using namespace arma;
using namespace std;

//----------------------------------------------------------
// Global Vars:
//----------------------------------------------------------
bool VERBOSE = false;
const int UPPER = 1;
const int LOWER = 0;

//----------------------------------------------------------
// Main Function
//----------------------------------------------------------
int main(int argc, char** argv)
{
    cout << endl << "-------------------------------------" << endl;
    cout << "LP Solver for General Simplex Problems" << endl;
    cout << "by Dave Coleman" << endl;
    cout << "-------------------------------------" << endl << endl;

	// Check if a parameter is being passed in - this is specifying what file to run
	if(argc > 1)
	{
		const char* filename = argv[1];
		cout << "Opening file: " << filename << endl << endl;
  	
		VERBOSE = true;
		dictionary s1 = readProblem(filename);
		
		s1 = solveLP(s1);
		
		outputResults(s1);
	}
	else // no filename was specified, just run all tests
	{
		cout << endl << "No filename passed as parameter to program. Running all tests" << endl;
		cout << endl << endl;
		runTests();
	}
	
	cout << endl;
	return EXIT_SUCCESS;
}
//----------------------------------------------------------
// Unit Testing, kinda
//----------------------------------------------------------
void runTests()
{
	double total_time = 0, start, end; // for benchmarking
	int num_tests = 100000;

	dictionary s1;
	
	// Solve problem 1 -------------------------------------------------------------
	VERBOSE = false;

	

	cout << "\033[1;31m" << "*********************************************************************************";
	cout << endl << "Running Example 1 - from Chapter 9 of Vanderbei" << endl << "\033[0m";
	
    s1 = readProblem("tests/example1.txt");

	// Benchmark the solver:
	start = get_time();
	for(int i = 0; i < num_tests; ++i)
	{
		s1 = solveLP(s1);
		checkObjective(s1, 3); // Check obj value
	}
	end = get_time();
	
	// Show results:	
	cout << "Ran " << num_tests << " in " << (end - start) << " seconds. " << endl;
	total_time += end - start;
	outputResults(s1);	

	

	
	// Solve problem 2 -------------------------------------------------------------
	VERBOSE = false;
	cout << "\033[1;31m" << "*********************************************************************************";
	cout << endl << "Running Example 2 - from Prof General Intialization PDF" << endl << "\033[0m";
	
    s1 = readProblem("tests/example2.txt");

	// Benchmark the solver:
	start = get_time();
	for(int i = 0; i < num_tests; ++i)
	{
		s1 = solveLP(s1);
		checkObjective(s1, 118); // Check obj value
	}
	end = get_time();
	
	// Show results:	
	cout << "Ran " << num_tests << " in " << (end - start) << " seconds. " << endl;
	total_time += end - start;
	outputResults(s1);



	
	

	// Solve problem 3 -------------------------------------------------------------
	VERBOSE = false;
	cout << "\033[1;31m" << "*********************************************************************************";
	cout << endl << "Running Example 3 - from Prof's Matlab Code, Example 1" << endl << "\033[0m";

    s1 = readProblem("tests/example3.txt");
	// Benchmark the solver:
	start = get_time();
	for(int i = 0; i < num_tests; ++i)
	{
		s1 = solveLP(s1);
		checkObjective(s1, 10.6667); // Check obj value
	}
	end = get_time();
	
	// Show results:	
	cout << "Ran " << num_tests << " in " << (end - start) << " seconds. " << endl;
	total_time += end - start;
	outputResults(s1);




	
	// Solve problem 4 -------------------------------------------------------------
	VERBOSE = false;
	cout << "\033[1;31m" << "*********************************************************************************";
	cout << endl << "Running Example 4 - from Prof's Notes Lecture 8, Slide 25" << endl << "\033[0m";

    s1 = readProblem("tests/example4.txt");

	// Benchmark the solver:
	start = get_time();
	for(int i = 0; i < num_tests; ++i)
	{
		s1 = solveLP(s1);
		checkObjective(s1, 3); // Check obj value
	}
	end = get_time();
	
	// Show results:	
	cout << "Ran " << num_tests << " in " << (end - start) << " seconds. " << endl;
	total_time += end - start;
	outputResults(s1);



	

	// Solve problem 5 -------------------------------------------------------------
	VERBOSE = false;
	cout << "\033[1;31m" << "*********************************************************************************";
	cout << endl << "Running Example 5 - from Chvatal Chapter 3 page 39" << endl << "\033[0m";
	
    s1 = readProblem("tests/example5.txt");


	// Benchmark the solver:
	start = get_time();
	for(int i = 0; i < num_tests; ++i)
	{
		s1 = solveLP(s1);
		checkObjective(s1, 4); // Check obj value
	}
	end = get_time();
	
	// Show results:	
	cout << "Ran " << num_tests << " in " << (end - start) << " seconds. " << endl;
	total_time += end - start;
	outputResults(s1);
	
	
	cout << endl << (5*num_tests) << " tests completed in " << total_time << " seconds." << endl;
}
//-------------------------------------------------------------------------------------------
// Check that the dictionary result is correct
//-------------------------------------------------------------------------------------------
void checkObjective(dictionary s1, double answer)
{
	// Check obj value
	if(fabs(s1.objvalue - answer) > .001)
	{
		cout << endl << "Incorrect answer found!" << endl;
		cout << "Solver Answer: " << s1.objvalue << endl;
		cout << "Correct Answer: " << answer << endl;
		throw;
	}
}
//-------------------------------------------------------------------------------------------
// Solve Linear Program
//-------------------------------------------------------------------------------------------
dictionary solveLP(dictionary s1)
{
	s1 = setup(s1); // setup variable tracking, etc
	
	// Output the Initial Dictionary
	if(VERBOSE)
	{
		cout << "Intial Dictionary" << endl;	
		printDictionary(s1);
	}
	
    // Check if its optimal
    if(isOptimal(s1))
	{
		if(VERBOSE)
			cout << "This LP problem is optimal!" << endl;
		
		return s1;
	}
	else
	{
		if(VERBOSE)
			cout << "Not optimal." << endl;
	}

	if(isFeasible(s1))
	{
		if(VERBOSE)
			cout << "This LP problem is feasible" << endl << endl;
	}
	else
	{
		if(VERBOSE)
			cout << "Not feasible." << endl << "Beginning initialization phase." << endl;

		s1 = initialize(s1);

		// Check if feasible now, just to be sure :-)
		if( !isFeasible(s1))
		{
			cout << "Still not feasible. Fail." << endl;
			throw;
		}

		// Check if optiaml, then we are done
		if(isOptimal(s1))
		{
			cout << "This LP problem is optimal!" << endl;
			return s1;
		}		
	}

	// Run simplex and return the result
	return simplex(s1);
}

dictionary simplex(dictionary s1)
{
	// Start Simplex -------------------------------------------------
	int step = 1; // keep track of how many pivots we do
	bool useBland = false;
	
	while(true)
	{
		// Pivot
		s1 = pivot(s1, useBland);

		// Check if dictionary is optimal
		if(isOptimal(s1))
		{
			cout << "Optimal Dictionary" << endl;
			if(VERBOSE)
				printDictionary(s1);
			break;
		}
		// Ensure termination by running Bland's rule after 50th step
		else if(step == 50)
		{
			useBland = true;
		}
		// Check for cycling or just bugs
		else if(step > 1000)
		{
			cout << "Possible cycling occuring, its been 1000 iterations. Ending" << endl;
			if(VERBOSE)
				printDictionary(s1);
			break;
		}
		else if(VERBOSE)		// Output the dictionary
		{
			cout << "Simplex Step " << step << endl;
			printDictionary(s1);
		}
		++step;
		
	}
	
    return s1;
}
//-------------------------------------------------------------------------------------------
// Setup remaining details of general form dictionary 
//-------------------------------------------------------------------------------------------
dictionary setup(dictionary s1)
{
	// Initialize the size of the nonbasic values matrix
    s1.nonbasic_values.set_size(1, s1.nonbasic.n_cols);
	
	// Decide how to rest the nonbasic variables
	for(int col = 0; col < int(s1.nonbasic.n_cols); ++ col)
	{
		// Rule 1: have bound rest on non-inifinte option
		if( !is_finite(s1.nonbasic_lower(0, col) ) ) // lower bound is negative inifinte
		{
			s1.nonbasic_values(0, col) = 1; // set to upper bound
		}
		else if( !is_finite(s1.nonbasic_upper(0, col) ) ) // upper bound is inifinte
		{
			s1.nonbasic_values(0, col) = 0; // set to lower bound
			// TODO: what if both lower and upper are infinite?
		}
		// Rule 2: Choose bound which maximizes the objective
		else
		{
			// so, if coefficient is negative choose lower, otherwise opposite
			if( s1.nonbasic(0, col) < 0 ) // is negative
			{
				s1.nonbasic_values(0, col) = 0; // set to lower
			}
			else
			{
				s1.nonbasic_values(0, col) = 1; // set to upper
			}
		}
	}

	// Setup the basic var values cache
	s1.basic_values.set_size(s1.basic.n_rows, 1);
	
	// Setup the variable tracking 
	// Nonbasic Vars:       0 to (n-1)
	// Basic Vars:          n to (n+m-1)
	// Auxillary Vars e:    (n+m) to (n+m+n_es-1)
    s1.nonbasic_vars.set_size(1, s1.nonbasic.n_cols);
    s1.basic_vars.set_size(s1.basic.n_rows, 1);
		
	for(int col = 0; col < int(s1.nonbasic.n_cols); ++col)
		s1.nonbasic_vars(0, col) = col; // index the rows
	
	for(int row = 0; row < int(s1.basic.n_rows); ++row)
		s1.basic_vars(row, 0) = row + int(s1.nonbasic.n_cols); // index the rows

	// Calculate the slack variables
	s1 = calculateSlack(s1, false);
	
	return s1;
}
//-------------------------------------------------------------------------------------------
// Get the row index of the corresponding leaving vaiable
//-------------------------------------------------------------------------------------------
void getLeavingVar(dictionary s1, int entering_var_index, int& leaving_var_index,
				   int& leaving_var_bound, bool useBland, bool& doFlip)
{
	// decide which constraint bounds it the most ( to the lowest value)
	// in other words, find t < a where a is the smallest
	double t = 0;
	int t_bound = -1;
	
	doFlip = false; // by default do not flip the bound
	
	// store the constraint that limits the obj val the most
	double smallest_const = numeric_limits<double>::infinity();
	int    smallest_const_index = 0;
	int    smallest_const_bound = -1; // 1 or 2 for lower or upper bound on smallest
	
	// loop through all constraints
	for(int row = 0; row < int(s1.basic.n_rows); ++row)
	{
		// Decide what constraint is active and lowest -------------------		

		//cout << "ROW " << row << " NONBASIC BOUNDED ON " << s1.nonbasic_values(0, entering_var_index)
		//	 << " COEFF " << s1.basic(row, entering_var_index) << endl;
		
		// CASE 1: Entering variable is on UPPER bound AND entering variable row coefficient is POSITIVE
		if( s1.nonbasic_values(0, entering_var_index) == UPPER &&
			s1.basic(row, entering_var_index) > 0)
		{
			// LOWER bound of basic row - current value of row
			t = s1.basic_values(row, 0) - s1.basic_lower(row, 0);
			t_bound = 0; // lower bound
			//cout << "CASE 1" << endl;
		}
   		// CASE 2: Entering variable is on LOWER bound AND entering variable row coefficient is POSITIVE
		else if( s1.nonbasic_values(0, entering_var_index) == LOWER &&
				 s1.basic(row, entering_var_index) > 0 )
		{
			// UPPER bound of basic row - current value of row  
			t = s1.basic_upper(row, 0) - s1.basic_values(row, 0);
			t_bound = 1; // upper bound			
			//cout << "CASE 2" << endl;			
		}
		// CASE 3: Entering variable is on LOWER bound AND entering variable row coefficient is NEGATIVE
		else if( s1.nonbasic_values(0, entering_var_index) == LOWER &&
			s1.basic(row, entering_var_index) < 0 )
		{
			// LOWER bound of basic row - current value of row
			t = s1.basic_lower(row, 0) - s1.basic_values(row, 0);	
			t_bound = 0; // lower bound		
			//cout << "CASE 3" << endl;			
		}
		// CASE 4: Entering variable is on UPPER bound AND entering variable row coefficient is NEGATIVE
		else if( s1.nonbasic_values(0, entering_var_index) == UPPER &&
			s1.basic(row, entering_var_index) < 0 )
		{
			// UPPER bound of basic row - current value of row 
			t = s1.basic_values(row, 0) - s1.basic_upper(row, 0);
			t_bound = 1; // upper bound						
		}
		// CASE 5: entering variable row coefficient is zero
		else if( s1.basic(row, entering_var_index) == 0)
		{
			// Not eligible to leave
			//cout << "CASE 5" << endl;
			continue; // move to next row
		}
		else
   		{
			cout << "No case found!" << endl;
			throw; //this shouldn't happen
		}
		
		// Divide by coefficient of current row
		t = t /  s1.basic(row, entering_var_index);
		
		//cout << endl << "Min Contraint: t <= " << t << " on row " << row << endl << endl;

		// now decide if this is the smallest value, or just the first one if bland is in effect
		if( t < smallest_const || useBland )
		{
			smallest_const = t;
			smallest_const_index = row;

			// The constraint is the same one used to calculate t
			// Except check to make sure said bound is not inifinity TODO: is this right?
			smallest_const_bound = t_bound;

			// Stop looking for other potential entering vars if bland's rule is in effect
			if(useBland)
			{
				break;
			}

			//cout << "smallest const = " << t << " index " << smallest_const_index
			//	 << " bounded at " << smallest_const_bound << endl;
		}

	}


   	// Now check if a flip variable is possibly better
	// 1. Check if coeff of nonbasic var is NEGATIVE and is on UPPER bound OR
	// 2. Check if coeff of nonbasic var is POSITIVE and is on LOWER bound
	if( ( s1.nonbasic(0, entering_var_index) < 0 && s1.nonbasic_values(0, entering_var_index) == UPPER ) ||
     	( s1.nonbasic(0, entering_var_index) > 0 && s1.nonbasic_values(0, entering_var_index) == LOWER ) )
	{
		double flip_amount = s1.nonbasic_upper(0, entering_var_index)
			- s1.nonbasic_lower(0, entering_var_index);

		if( !is_finite(flip_amount) )
		{
			//cout << "Flip ignored because is inifinity = " << flip_amount << endl;
		}
		else
		{
			if(VERBOSE)
				cout << "Flip is possible, amount = " << flip_amount << "." << endl;

			if( flip_amount < smallest_const )
			{
				if(VERBOSE)				
					cout << "Flipping. Is least bound on increase of entering variable" << endl << endl;
				doFlip = true;
			}
			else
			{
				if(VERBOSE)				
					cout << "Flip lost to " << smallest_const << endl;

			}
		}
	}
	
	// check that a leaving variable was found
	if( smallest_const == numeric_limits<double>::infinity() )
	{
		cout << "No leaving variable found. Problem is unbounded." << endl;
		throw;
	}

	assert(smallest_const_bound > -1); // make sure this value was updated
	
	leaving_var_index = smallest_const_index;
	leaving_var_bound = smallest_const_bound;
}

//-------------------------------------------------------------------------------------------
// Get the row index of the next entering variable
//-------------------------------------------------------------------------------------------
void getEnteringVar(dictionary s1, bool useBland, int& entering_var_index)
{
	// NEW METHOD:
	// c*x is a candidate for leaving if:
	//  1) x is on its lower bound AND c  > 0
	//  OR
	//  2) x is on its upper bound AND c < 0
	//
	// Choose the candidate that increases z the most
	//
	// Also look for flip variables

	double largest_coef = -1*numeric_limits<double>::infinity(); // keep track of the largest found coefficient
	int largest_coef_index = -1; // init with a not found flag	

	for(int col = 0; col < int(s1.nonbasic.n_cols); ++col)
	{
		// 1) Check if coefficeint is positive and at lower bound
		if(s1.nonbasic(0, col) > 0 && s1.nonbasic_values(0, col) == LOWER)
		{
			// is candidate to leave
			if(fabs(s1.nonbasic(0, col)) > largest_coef)
			{
				largest_coef = fabs(s1.nonbasic(0, col));
				largest_coef_index = col;

				// Stop looking for other potential entering vars if bland's rule is in effect
				if(useBland)
				{
					break;
				}				
			}
		}
		// 2) Check if coefficient is negative and at upper bound
		else if(s1.nonbasic(0, col) < 0 && s1.nonbasic_values(0, col) == UPPER)
		{
			
			// is candidate to leave
			if(fabs(s1.nonbasic(0, col)) > largest_coef)
			{
				largest_coef = fabs(s1.nonbasic(0, col));
				largest_coef_index = col;

				// Stop looking for other potential entering vars if bland's rule is in effect
				if(useBland)
				{
					break;
				}								
			}			
		}
	}

	// Check if no entering variable found. This should not happen
	if(largest_coef == -1*numeric_limits<double>::infinity())
	{
		cout << "No entering variable found. This function is not feasible." << endl;
		throw;
	}
	
	// Return results by reference
	entering_var_index = largest_coef_index;
}
//-------------------------------------------------------------------------------------------
// Check if current general form is optimal
//-------------------------------------------------------------------------------------------
bool isOptimal(dictionary s1)
{
	// NEW METHOD:
	// Check that all coefficients in the obj function are:
	//  1) negative and at lower bound
	//  OR
	//  2) positive and at upper bound
	//
	// If this isnt true for all coefficients, prob is not optimal
	
	for(int col = 0; col < int(s1.nonbasic.n_cols); ++col)
	{
		// 1) Check if coefficeint is negative and at lower bound
		if(s1.nonbasic(0, col) < 0 && s1.nonbasic_values(0, col) == LOWER)
		{
			// this var is as optimal as it can get
		}
		// 2) Check if coefficient is positive and at upper bound
		else if(s1.nonbasic(0, col) > 0 && s1.nonbasic_values(0, col) == UPPER)
		{
			// this var is as optimal as it can get
		}
		else if(s1.nonbasic(0, col) == 0)
		{
			// 0 means optimal
		}
		else
		{
			return false;
		}

	}
	
	// OLD:
	// Check if all coefficients in objective function are negative
	// if not negative, check if the positive coefficient is at its upper bound
	// if not at upper bound, the solution is not optimal

		
		/* Check if coefficient is positive
		if(s1.nonbasic(0, col) >= 0) 
		{
			// is positive
			
			// Check if positive coefficient is at its upper bound
			if( s1.nonbasic_values(0, col) != 1 )
			{
				// is NOT at upper bound
				return false;
			}

		}
		*/

	
	return true;
}
//-------------------------------------------------------------------------------------------
// Calculates current value of each row
// checkObjValue is a flag that should only be false if we are doing auxillary problem
//-------------------------------------------------------------------------------------------
dictionary calculateSlack(dictionary s1, bool checkObjValue)
{
	//Calculate values of slack variables, Wn...
	
	double value; // holds the current row's solution value
	
	// Loop through each constraint
	for(int row = 0; row < int(s1.basic.n_rows); ++row)
	{
		value = 0; //reset the row value amount
		
	    for( int col = 0; col < int(s1.basic.n_cols); ++col)
		{
			// TODO: turn this into matrix math stuff.... ?
			value += s1.basic(row, col) * getNonbasicVal(s1, col);			
		}
		s1.basic_values(row, 0) = value;
	}

	value = 0; // reuse this var
	
	// Calculate the objective value
	for(int col = 0; col < int(s1.nonbasic.n_cols); ++col)
	{
		value = value + s1.nonbasic(0, col) * getNonbasicVal(s1, col);
	}

	// Check if degenerate
	if(VERBOSE && s1.objvalue == value)
	{
		cout << "The objective value did not change. Degeneracy." << endl << endl;
	}

	// Check for error
	if( s1.objvalue > value && checkObjValue)
	{
		cout << endl << "Old ObjVal = " << s1.objvalue << endl;
		cout << "New ObjVal = " << value << endl;
		cout << "OBJECTIVE VALUE DECREASED. SOMETHING WENT WRONG" << endl;
		throw;
	}

	//Update value
	s1.objvalue = value;
	
	return s1;
}
//-------------------------------------------------------------------------------------------
// Check if dictionaty is feasible
//-------------------------------------------------------------------------------------------
bool isFeasible(dictionary s1)
{
	// Loop through each constraint
	for(int row = 0; row < int(s1.basic.n_rows); ++row)
	{
		// Now check if value is within constraint bounds
		if( ! ( s1.basic_values(row, 0) >= s1.basic_lower(row, 0) &&
			    s1.basic_values(row, 0) <= s1.basic_upper(row, 0) ) )
		{
			return false; // constraint not satisfied
		}
	}

	// Is within bounds
	return true;
}
//-------------------------------------------------------------------------------------------
// Do the actual pivot
//-------------------------------------------------------------------------------------------
dictionary pivot(dictionary s1, bool useBland)	
{
	// Step 0: Find Entering and leaving variables ---------------------------------------
	
    // Choose the entering variable
	int entering_var_index;
	
	getEnteringVar(s1, useBland, entering_var_index);
	if(VERBOSE)
		cout << "Entering variable: " << resolveVarName(s1, s1.nonbasic_vars[entering_var_index]) << endl;

	// Choose the leaving variable
	int leaving_var_index, leaving_var_bound;
	bool doFlip = false; // passed as reference variable and updated if pivot is not best step
	getLeavingVar(s1, entering_var_index, leaving_var_index, leaving_var_bound, useBland, doFlip);

	// Check if a flip is better than the leaving variable
	if( doFlip )
	{
		// Forget the pivot, we are flipping the bound of the entering variable and calling it a day
		s1.nonbasic_values(0, entering_var_index) = !s1.nonbasic_values(0, entering_var_index);
   		s1 = calculateSlack(s1, true);
		
		return s1;
	}
	
	if(VERBOSE)
		cout << "Leaving variable: " << resolveVarName(s1, s1.basic_vars[leaving_var_index]) << endl << endl;
	
	// Step 1: Create Replacement Rule --------------------------------------------------
	
	// First create the replacement rule
	mat replacement_rule =s1.basic.row(leaving_var_index);

	// Replace the entering var coeff with a negative 1
	replacement_rule(0, entering_var_index) = -1; // b/c we've subtracted from other side
		
	// Divide whole row by entering var coefficient
	replacement_rule = replacement_rule / (-1 * s1.basic(leaving_var_index, entering_var_index));

	// Store the replacement rule back into the original row
	s1.basic.row(leaving_var_index) = replacement_rule;

	// Swap the variable names
	int new_basic_var = s1.nonbasic_vars(0, entering_var_index);
	s1.nonbasic_vars(0, entering_var_index) = s1.basic_vars(leaving_var_index, 0);
	s1.basic_vars(leaving_var_index, 0) = new_basic_var;

	// Step 2: Substitute into Obj Function ------------------------------------------
	
	// Create replacement rule for obj function (step 1 of substiting in rep rule)
	mat rep_rule1 = replacement_rule * s1.nonbasic(0, entering_var_index);

	// Change entering variable position in nonbasic function to 0, temporarily
	s1.nonbasic(0, entering_var_index) = 0;

	// Now add the modified nonbasic rule with the custom replacement rule
	s1.nonbasic = s1.nonbasic + rep_rule1;

	// Swap upper and lower bounds from basic to nonbasic
	double new_nonbasic_upper = s1.basic_upper(leaving_var_index, 0);
	double new_nonbasic_lower = s1.basic_lower(leaving_var_index, 0);
	s1.basic_upper(leaving_var_index, 0) = s1.nonbasic_upper(0, entering_var_index);
	s1.basic_lower(leaving_var_index, 0) = s1.nonbasic_lower(0, entering_var_index);	
	s1.nonbasic_upper(0, entering_var_index) = new_nonbasic_upper;
	s1.nonbasic_lower(0, entering_var_index) = new_nonbasic_lower;

	// Step 3: Change the constraints the non-basic variables are resting on ----------

	// If the leaving_var_bound was set to lower bound, move the nonbasic var bound to lower
	if(leaving_var_bound == LOWER) // lower
	{
		s1.nonbasic_values(0, entering_var_index) = 0; // switch to lower
	}
	else
	{
		s1.nonbasic_values(0, entering_var_index) = 1; // switch to upper
	}
	
	// Step 4: Substitute into rest of constraints ------------------------------------

	// Loop through every row
	for(int row = 0; row < int(s1.basic.n_rows); ++row)
	{
		// Check if this row is not the replacement rule row (leaving row)
		if(row != leaving_var_index)
		{
			//cout << "Appling rep rule to row" << row  << endl;
			// It is not, do substitution
			rep_rule1 = replacement_rule * s1.basic(row, entering_var_index);
			
			// Replace entering var location with 0
			s1.basic(row, entering_var_index) = 0;

			// Now add the custom replacement rule to the row
			s1.basic.row(row) = s1.basic.row(row) + rep_rule1;
			
		}
	}

	// Step 5: Update slack variable amount/row solutions
	s1 = calculateSlack(s1, true); 

	return s1;
}
//-------------------------------------------------------------------------------------------
// Intialize the problem through dualization
//-------------------------------------------------------------------------------------------
dictionary initialize(dictionary s1)
{
	// Create a copy of the s1 general form dictionary
	// and perform initialization until problem becomes feasible

	dictionary s2 = s1;
	
	// Introduce variables e1 to ek to force infeasible variables back into appropriate ranges

	double value; // holds the current row's solution value. this tells us the sign for var e

	// Assume we will be adding auxillary variables, so go ahead and clear out objective function
	s2.nonbasic.fill(0); // set all to zero

	// Loop through each constraint
	for(int row = 0; row < int(s2.basic.n_rows); ++row)
	{
		value = 0; //reset the row value amount
		
	    for( int col = 0; col < int(s2.basic.n_cols); ++col)
		{
			value += s2.basic(row, col) * getNonbasicVal(s2, col);
		}

		// Now check if value is within constraint bounds
		if( ! ( value >= s2.basic_lower(row, 0) &&
			    value <= s2.basic_upper(row, 0) ) )
		{
			// constraint not satisfied. so now we add e variable
			if(VERBOSE)
				cout << "Out of bounds constraint found on row " << row
					 << " with value " << value << "." << endl << endl;

			// Add column at right of matrix with all zeros
			s2.basic.insert_cols(s2.basic.n_cols, 1);
			
			// Add new column to all nonbasic matricies
			s2.nonbasic.insert_cols(s2.nonbasic.n_cols, 1);
			s2.nonbasic_lower.insert_cols(s2.nonbasic_lower.n_cols, 1);
			s2.nonbasic_upper.insert_cols(s2.nonbasic_upper.n_cols, 1);
			s2.nonbasic_values.insert_cols(s2.nonbasic_values.n_cols, 1);
			s2.nonbasic_vars.insert_cols(s2.nonbasic_vars.n_cols, 1);

			// Set new objective to negative of auxillary variable
			s2.nonbasic(0, s2.nonbasic.n_cols - 1) = -1;

			// Always set lower bound of ex to 0
			s2.nonbasic_lower(0, s2.nonbasic_lower.n_cols - 1) = 0;

			// Set nonbasic constraint to rest on upper bound
			s2.nonbasic_values(0, s2.nonbasic_values.n_cols - 1) = 1;

			// Set nonbasic var names
			s2.nonbasic_vars(0, s2.nonbasic_vars.n_cols - 1) = s2.nonbasic_vars.n_cols + s2.basic_vars.n_rows - 1;

			// Check if TOO SMALL
			if( value < s2.basic_lower(row, 0) )
			{
				s2.basic(row, s2.basic.n_cols - 1) = 1; // ADD ek to row

				// set upper bound equal to LOWER LIMIT - SOLUTION
				s2.nonbasic_upper(0, s2.nonbasic_upper.n_cols - 1) =
					s2.basic_lower(row, 0) - value;
			}

			// CHECK if TOO LARGE
			if( value > s2.basic_upper(row, 0) )
			{
				s2.basic(row, s2.basic.n_cols - 1) = -1; // SUBTRACT ek from row

				// set upper bound equal to SOLUTION - UPPER LIMIT
				s2.nonbasic_upper(0, s2.nonbasic_upper.n_cols - 1) =
					value - s2.basic_upper(row, 0);
			}

		}
	}

	// Update slack variable amount/row solutions
	s2 = calculateSlack(s2, false);		
	
	// We now have a dictionary ready for the initialization phase:
	if(VERBOSE)
	{
		cout << "New Auxillary Problem " << endl;		
		printDictionary(s2);
	}
	
	s2 = simplex(s2);

	// We should now have an optimal intial dictionary
	// Now we check that the original dictionary is feasible
	// Check that aux LP obj val is 0
	if( s2.objvalue != 0)
	{
		// Orig LP Problem is not feasible
		cout << "New dictionary is not feasible!" << endl;
		throw;
	}
	
	// Next we remove the nonbasic independent variables

	// Calculate the auxillary id threshold:
	int aux_id_min = s2.n + s2.m;
	
	for(int col = 0; col < int(s2.nonbasic.n_cols); ++col)
	{
		// Check if var id is within aux range
		if( s2.nonbasic_vars(0, col) >= aux_id_min )
		{
			// this is an aux variable

			// remove column from all matricies/vectors
			if(VERBOSE)
				cout << "Removing column " << col << "." << endl << endl;
			
			s2.nonbasic.shed_col(col);
			s2.nonbasic_values.shed_col(col);
			s2.nonbasic_lower.shed_col(col);
			s2.nonbasic_upper.shed_col(col);
			s2.nonbasic_vars.shed_col(col);
			s2.basic.shed_col(col);

		}
	}
	
	// if any of aux vars are basic, it cannot be removed
	// instead we set both its bounds to 0
	for(int row =0; row < int(s2.basic.n_rows); ++row)
	{
		// Check if var id is within aux range
		if( s2.basic_vars(row, 0) >= aux_id_min )
		{
			// this is an aux variable

			// set bounds to 0
			s2.basic_lower(row, 0) = 0;
			s2.basic_upper(row, 0) = 0;			
		}
	}
	
	// re-add original objective function
	s1 = combineObjFunc(s1, s2);

	// Update all the row totals and obj value
	s1 = calculateSlack(s1, false);
	
	if(VERBOSE)
	{
		cout << "After Converting Dictionary Back To Original LP" << endl;
		printDictionary(s1);

		cout << endl << "Done with initilization." << endl;		
	}
	
	return s1;
}
//-------------------------------------------------------------------------------------------
// Combine original objective function to auxilary dictionary
//-------------------------------------------------------------------------------------------
dictionary combineObjFunc(dictionary s1, dictionary s2)
{
	// Move obj function from s1 into s2. Return s2
	if(VERBOSE)
		cout << "COMBING OBJECTIVE FUNCTIONS " << endl << endl << endl << endl;
	
	bool haveAdded; // keep track of each var status
	
	// First reset s2 obj function to zero
	s2.nonbasic.fill(0);

	// Now loop through every col of old obj function and add components to new obj function
	for(int col = 0; col < int(s1.nonbasic.n_cols); ++col)
	{
		haveAdded = false;

		// check if this variable is a nonbasic var in the aux dictionary
		// if it is, simply add it
		for( int col2 = 0; col2 < int(s2.nonbasic.n_cols); ++col2)
		{
			if( s1.nonbasic_vars(0, col) == s2.nonbasic_vars(0, col2) )
			{
				
				// this nonbasic variable appears non basic in both dictionarys, direct add
				s2.nonbasic(0, col2) = s2.nonbasic(0, col2) + s1.nonbasic(0, col);
				haveAdded = true;

				break;
			}
		}

		if( !haveAdded )
		{
			// this nonbasic variable does not appear in both dictionaries
			// multiply instead the coefficient of original obj function times the basic row
			// add then add this basic row to the new objective function

			// First, find the basic row that has the same variable id. This is so stupid complicated.
			for( int row = 0; row < int(s2.basic.n_rows); ++row)
			{
				if( s2.basic_vars(row, 0) == s1.nonbasic_vars(0, col) )
				{
					// the original obj function var has been matched to a nonbasic row in the aux LP

					mat tempRow = s2.basic.row(row);

					// Multiply this row by the original obj variable's coefficient
					tempRow = tempRow * s1.nonbasic(0, col);

					// Add this temp row to the new objective function
					s2.nonbasic = s2.nonbasic + tempRow;

				}

			}
			
		}

	}
	
	return s2;
}
//-------------------------------------------------------------------------------------------
// Compare 2 Dictionaries
//-------------------------------------------------------------------------------------------
bool dictionaryIsEqual(dictionary d1, dictionary d2)
{
	// Compare nonbasic matrix
	umat compareMatrix = (d1.nonbasic == d2.nonbasic);
	
	if( accu(compareMatrix) == double( d1.nonbasic.n_cols * d1.nonbasic.n_rows ) )
	{
		// Compare basic matrix
		compareMatrix = (d1.basic == d2.basic);
	
		if( accu(compareMatrix) == double( d1.basic.n_cols * d1.basic.n_rows ) )
		{
			if(VERBOSE)
			{
				cout << "Answer is correct" << endl;
			}
			return true;
		}
		else
		{
			cout << "Basic matrix not the same." << endl;
			d1.basic.print("D1 Basic");
			d2.basic.print("D2 Basic");
			return false;
		}
	}
	else
	{
		cout << "Non basic matrix not the same." << endl;
		d1.nonbasic.print("D1 Nonbasic");
		d2.nonbasic.print("D2 Nonbasic");
		
		return false;
	}
}
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
