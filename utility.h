#include <armadillo>

using namespace arma;
using namespace std;

/*
  Objective function vector c⃗ 
  The matrix A
  The row bounds: a⃗ ,b⃗  (some of the entries can be Inf or -Inf ).
  The variable bounds: l⃗ ,u⃗ .
  The input format is quite simple text format.

  m,n %% #Constraints #Variables
  c1, ...,cn %% objective coeffs
  a11, ..., a1n %% matrix first row
  ...
  am1, ..., amn %% matrix last row
  a1 ... am %% Lower bound vector: can have -Inf entries
  b1 ... bm %% Upper bound vector: can have Inf entries
  l1... ln %% Var lower bounds can have -Inf entries
  u1... un %% Var upper bounds can have Inf entries
*/
struct genForm
{
	int m; // number constraints
	int n; // number variables
	// Objective function
	mat c;
	// Matrix A
	mat A; // a00 to a0n, am0 to amn
	// Row bounds:
	mat a, b;
	// Variable bounds:
	mat l, u;
};

struct simplex
{
	mat nonbasic; // the row of the non basic variables in the objective function
	mat nonbasic_values; // bool row -   0 = lower and     1 = upper
	mat nonbasic_lower; // row of lower bounds for nonbasic vars
	mat nonbasic_upper; // row of upper bounds for nonbasic vars
	mat nonbasic_vars; // tracks the names of the variables as they move with simplex
	
	mat basic; // matrix of basic variables that are the constraints
	mat basic_lower; // col of lower bounds for basic vars
	mat basic_upper; // col of upper bounds for basic vars
	mat basic_vars; // tracks the names of the variables as they move with simplex
};

// Read a problem from file
genForm readProblem();

// Read a row and convert to matrix
arma::mat readRow(int mn, std::string line);

// Convert a cell
double convertCell(std::string &s);

// Output to screen a gen form equation
void printGenForm(genForm genForm1);

// Output to screen a formatted simplex step
void printSimplexStep(simplex simplex1);

// format positive numbers less than 10
string tabber(double num);
