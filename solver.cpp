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

using namespace arma;
using namespace std;

//----------------------------------------------------------
// Global Vars:
//----------------------------------------------------------


//----------------------------------------------------------
// Main Function
//----------------------------------------------------------
int main(int argc, char** argv)
{
    cout << endl;
    cout << "-------------------------------------" << endl;
    cout << "LP Solver for General Simplex Problems" << endl;
    cout << "by Dave Coleman" << endl;
    cout << "-------------------------------------" << endl << endl;

    genForm genForm1 = readProblem();

    // Output structure:
    //printGenForm(genForm1);

    // TODO: check for feasibility - if all wi's are between upper and lower bounds

    // Convert genForm1 to simplex form
    simplex simplex1;

    // Nonbasic/decision variables:
    simplex1.nonbasic = genForm1.c;
    simplex1.nonbasic_lower = genForm1.l;
    simplex1.nonbasic_upper = genForm1.u; 
    simplex1.nonbasic_values.set_size(1, genForm1.n);
    simplex1.nonbasic_values.fill(0.0);
	// TODO: smarter way of choosing upper or lower for init nonbasic values?

    // Basic/Dependent Variables:
    simplex1.basic = genForm1.A;
    simplex1.basic_lower = trans(genForm1.a);
    simplex1.basic_upper = trans(genForm1.b);

	// Setup the variable tracking
	// the nonbasic vars are numbered 0 to (n-1) and the basic vars
	// are numbered n to (n+m-1)
	for(int row = 0; row < int(simplex1.nonbasic.n_rows); ++row)
		simplex1.nonbasic_vars(row, 0) = row; // index the rows

	for(int col = 0; col < int(simplex1.basic.n_cols); ++col)
		simplex1.nonbasic_vars(col, 0) = col + simplex1.nonbasic.n_rows; // index the cols

	// Output the dictionary
    cout << "Simplex Step 0" << endl;
    printSimplexStep(simplex1);

    // Check if its optimal
    if(isOptimal(simplex1))
	{
		cout << "This LP problem is optimal!" << endl;
	}
	else
	{
		cout << "Not optimal" << endl;
	}

	if(isFeasible(simplex1))
	{
		cout << "This LP problem is feasible" << endl;
	}
	else
	{
		cout << "Not feasible" << endl;
	}

	// Choose the entering variable
	int entering_var_index = getEnteringVar(simplex1);
	cout << "Entering variable: " << entering_var_index << endl;
	// Choose the leaving variable
	int leaving_var_index = getLeavingVar(simplex1, entering_var_index);
	cout << "Leaving variable: " << leaving_var_index << endl;

	// Pivot
	simplex1 = pivot(simplex1, entering_var_index, leaving_var_index);
	
	cout << endl;
    return EXIT_SUCCESS;
}
//-------------------------------------------------------------------------------------------
// Get the row index of the corresponding leaving vaiable
//-------------------------------------------------------------------------------------------
int getLeavingVar(simplex s1, int entering_var_index)
{
	// decide which constraint bounds it the most ( to the lowest value)
	// in other words, find t < a where a is the smallest

	double lower, upper, coeff = 0, other_coeff_total = 0, real_constraint;
	
	// store the constraint that limits the obj val the most
	double smallest_const = numeric_limits<double>::infinity();
	int    smallest_const_index = 0;

	// loop through all constraints
	for(int row = 0; row < int(s1.basic.n_rows); ++row)
	{
		// unload the data just so i can think!
		lower = s1.basic_lower(row, 0);
		upper = s1.basic_upper(row, 0);

		// loop through all coefficents in constraint
		for(int col = 0; col < int(s1.basic.n_cols); ++col)
		{
			// check if current constraint is the entering var
			if(col == entering_var_index)
			{
				// this is the entering variable coeff for this constraint
				coeff = s1.basic(row,col);
			}
			else
			{
				// calculate value of rest of constraint
			    other_coeff_total += s1.basic(row,col) * getNonbasicVal(s1, col);
			}

		}
		
		// now decide which side of the constraint we care about
		// this is determined by the sign of the entering variable in this contraint
		if(coeff < 0) // use lower bound
		{
			// move other constraint stuff to the other side of equality sign
			real_constraint = lower - other_coeff_total;
		}
		else // use upper bound
		{
			real_constraint = upper - other_coeff_total;
		}
			
		// move entering var coefficient to the other side
		real_constraint = real_constraint / coeff;

		// now decide if this is the smallest value
		if( real_constraint < smallest_const )
		{
			smallest_const = real_constraint;
			smallest_const_index = row;
		}

	}
	// check that a leaving variable was found
	if( smallest_const == numeric_limits<double>::infinity() )
	{
		cout << "No leaving variable found. Problem is unbounded." << endl;
		throw;
	}
	
	return smallest_const_index;
}
//-------------------------------------------------------------------------------------------
// Get the row index of the next entering variable
//-------------------------------------------------------------------------------------------
int getEnteringVar(simplex s1)
{
	// Does the same thing as isOptimal, except it checks all the vars and picks the best
	// one to be the entering variable

    u32 largest_coef_index;
	double largest_coef = s1.nonbasic.max(largest_coef_index);

	//cout << "Largest Coef: " << largest_coef << endl;
	//cout << "Index: " << largest_coef_index << endl;

	// Check if not entering variable found. This should not happen
	if(largest_coef <= 0)
	{
		cout << "No entering variable found. This function is not feasible." << endl;
		throw;
	}
	
	return largest_coef_index;
}
//-------------------------------------------------------------------------------------------
// Check if current general form is optimal
//-------------------------------------------------------------------------------------------
bool isOptimal(simplex s1)
{
	// Todo: make sure this is the correct way to check for optimality
   
	// Check if all coefficients in objective function are negative
	// if not negative, check if the positive coefficient is at its upper bound
	// if not at upper bound, the solution is not optimal
	for(int col = 0; col < int(s1.nonbasic.n_cols); ++col)
	{
		// Check if coefficient is positive
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
	}
	

	return true;
}
//-------------------------------------------------------------------------------------------
// Check if dictionaty is feasible
//-------------------------------------------------------------------------------------------
bool isFeasible(simplex s1)
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

		cout << "Constraint Value = " << value << endl;
		
		// Now check if value is within constraint bounds
		if( ! ( value >= s1.basic_lower(row, 0) &&
			    value <= s1.basic_upper(row, 0) ) )
		{
			return false; // constraint not satisfied
		}
	}

	// Is within bounds
	return true;
}
//-------------------------------------------------------------------------------------------
// Return the bound that the nonbasic variable is currently resting on, or "active'
//-------------------------------------------------------------------------------------------
double getNonbasicVal(simplex s1, int index)
{
	if( s1.nonbasic_values(0, index) ) // value is on lower bound
		return s1.nonbasic_upper(0, index);
	else
		return s1.nonbasic_lower(0, index);
}


//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
