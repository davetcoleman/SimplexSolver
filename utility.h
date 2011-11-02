#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;




struct dictionary
{
	int n; // number of ORIGINAL nonbasic vars. used for printing variable names only
	int m; // number of ORIGINAL basic vars. used for printing variable names only

	mat nonbasic; // the row of the non basic variables in the objective function
	mat nonbasic_values; // bool row -   0 = lower and     1 = upper
	mat nonbasic_lower; // row of lower bounds for nonbasic vars
	mat nonbasic_upper; // row of upper bounds for nonbasic vars
	mat nonbasic_vars; // tracks the names of the variables as they move with simplex
	
	mat basic; // matrix of basic variables that are the constraints
	mat basic_lower; // col of lower bounds for basic vars
	mat basic_upper; // col of upper bounds for basic vars
	mat basic_vars; // tracks the names of the variables as they move with simplex
	mat basic_values; // caches the value of each row, mostly just for printing

	double objvalue;
};

// Read a problem from file
dictionary readProblem(const char* filename);

// Read a row and convert to matrix
arma::mat readRow(int mn, std::string line);

// Convert a cell
double convertCell(std::string &s);

// Output to screen a formatted simplex step
void printDictionary(dictionary s1);

// Convert a variable name index into a name string
string resolveVarName(dictionary s1, int var_index);

// format positive numbers less than 10
string tabber(double num);
