// LP Solver for General Simplex Problems
// CSCI Linear Programming - Programming Assignment 1
// Dave Coleman | david.t.coleman@colorado.edu
// 10/20/2011
//
// Utility functions

//----------------------------------------------------------
// Directives:
//----------------------------------------------------------
#include <fstream>
#include "utility.h"

using namespace arma;
using namespace std;

// trim from both ends
static inline std::string &trim(std::string &s);

// trim from start
static inline std::string &ltrim(std::string &s);

// trim from end
static inline std::string &rtrim(std::string &s);

//-------------------------------------------------------------------------------------------
/*
  Read the problem file in .txt format and convert to matrix form
   
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
//-------------------------------------------------------------------------------------------
dictionary readProblem(const char* filename)
{
	// Create new general form problem
	dictionary inputGenForm;

    // Read the actual file
    std::ifstream  data(filename);

    std::string line;

    // --------------------------------------------------------------------
    // The first line should have the number of constraints and variables
    // number constrains:
    std::getline(data,line);
    std::stringstream  lineStream(line);
    std::string        cell;

    // Get the number of constraints
    getline(lineStream,cell,',');
    inputGenForm.m = convertCell(cell);

    // Get the number of variables
    getline(lineStream,cell,',');
    inputGenForm.n = convertCell(cell);

    // --------------------------------------------------------------------
    // The second line has the objective coeffients

    std::getline(data,line);
    inputGenForm.nonbasic = readRow(inputGenForm.n, line);

    // --------------------------------------------------------------------
    // The next m lines have the matrix rows

    mat A(inputGenForm.m,inputGenForm.n);    // Set the size of matrix A

    for(int m = 0; m < inputGenForm.m; ++m)
    {
		getline(data,line);
		std::stringstream  lineStream(line);

		for(int n = 0; n < inputGenForm.n; ++n)
		{
			// Get next number
			getline(lineStream,cell,',');

			// Save next number to matrix
			A(m,n) = convertCell(cell);
		}
    }

    inputGenForm.basic = A;

    // --------------------------------------------------------------------
    // The fourth to last line is lower boud vector
    std::getline(data,line);
    inputGenForm.basic_lower = trans(readRow(inputGenForm.m, line));

    // --------------------------------------------------------------------
    // 3rd to last: uppper bound vector
    std::getline(data,line);
    inputGenForm.basic_upper = trans(readRow(inputGenForm.m, line));

    // --------------------------------------------------------------------
    // 2nd to last: var lower bound
    std::getline(data,line);
    inputGenForm.nonbasic_lower = readRow(inputGenForm.n, line);

    // --------------------------------------------------------------------
    // last: var upper bounds
    std::getline(data,line);
    inputGenForm.nonbasic_upper = readRow(inputGenForm.n, line);

    return inputGenForm;
}
//-------------------------------------------------------------------------------------------
// Read a row of width mn into a matrix
//-------------------------------------------------------------------------------------------
arma::mat readRow(int mn, std::string line)
{
      std::stringstream  lineStream(line);
      std::string cell;

      mat x(1,mn);

      for(int n = 0; n < mn; ++n)
      {
	// Get next number
	getline(lineStream,cell,',');

	// Save next number to matrix
	x(0,n) = convertCell(cell);
      }

      return x;
}

//-------------------------------------------------------------------------------------------
// Trim whitespace around string, convert to int or inifinity value
//-------------------------------------------------------------------------------------------
double convertCell(std::string &s) {
  s = trim(s);
  double value = atof(s.c_str());
  return value;
}
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
void printDictionary(dictionary s1)
{
    cout << "-----------------------------------------------------------------------------------" << endl;
	
	// Print variable locations
	cout << "\t \t |    ";
	for(int col = 0; col < int(s1.nonbasic_vars.n_cols); ++col)
	{
		cout << resolveVarName(s1, s1.nonbasic_vars(0, col)) << "  \t";
	}
	cout << endl;
	
    // Print lower bounds
    cout << "l \t \t |    ";
    for(int col = 0; col < int(s1.nonbasic_lower.n_cols); ++col)
    {
        // check if the nonbasic variable is currently on this bound
		if(s1.nonbasic_values(0,col) == 0)
		{
			cout << "\033[1;31m"<< tabber(s1.nonbasic_lower(0,col)) << "\033[0m \t";
			//cout << "["<< tabber(s1.nonbasic_lower(0,col)) << "] \t";			
		}
		else
		{
			cout << tabber(s1.nonbasic_lower(0,col)) << "\t";
		}
    }
    cout << endl;

    // Print upper bounds
    cout << "\t u \t |    ";
    for(int col = 0; col < int(s1.nonbasic_upper.n_cols); ++col)
    {
        // check if the nonbasic variable is currently on this bound
        if(s1.nonbasic_values(0,col) == 1)
		{
			cout << "\033[1;31m"<< tabber(s1.nonbasic_upper(0,col)) << "\033[0m \t";
			//cout << "["<< tabber(s1.nonbasic_upper(0,col)) << "] \t";			
		}
		else
		{
			cout << tabber(s1.nonbasic_upper(0,col)) << " \t";
		}
    }
    cout << endl;
    cout << "-----------------------------------------------------------------------------------" << endl;	

    
    // Print objective function
    cout << "\t\t |z=  ";
    for(int col = 0; col < int(s1.nonbasic.n_cols); ++col)
    {
      cout << tabber(s1.nonbasic(0,col)) << " \t";
    }
    cout << " =" << s1.objvalue << endl;

    cout << "-----------------------------------------------------------------------------------" << endl;

	// Print basic variables
    for(double row = 0; row < int(s1.basic.n_rows); ++row)
    {
        cout << s1.basic_lower(row,0) << " \t ";
        cout << s1.basic_upper(row,0) << " \t |";
		cout << resolveVarName(s1, s1.basic_vars(row, 0)) << "  ";
	 
		for(int col = 0; col < int(s1.basic.n_cols); ++col)
		{
			cout << tabber(s1.basic(row, col)) << " \t ";
		}

		// Print basic values
		cout << "=" << tabber(s1.basic_values(row, 0)) << endl;

    }
    cout << "-----------------------------------------------------------------------------------" << endl;	


}
//-------------------------------------------------------------------------------------------
// Print the x values at end of program
//-------------------------------------------------------------------------------------------
void outputResults(dictionary s1)
{
	cout << endl << "Final Results:" << endl;

	mat result(1, s1.n); // store the x values only
	int nonbasic_id_max = s1.n; // only store ids that are less than this number

	// Output nonbasic values
    for(int col = 0; col < int(s1.nonbasic.n_cols); ++col)
	{
		if( s1.nonbasic_vars(0, col) < nonbasic_id_max )
		{
			result(0, s1.nonbasic_vars(0, col)) = getNonbasicVal(s1, col);
		}
		//cout << resolveVarName(s1, s1.nonbasic_vars(0, col)) << " = " << getNonbasicVal(s1, col) << endl;
	}

   	// Output basic values
	for(int row = 0; row < int(s1.basic.n_rows); ++row)
	{
		if( s1.basic_vars(row, 0) < nonbasic_id_max )
		{
			result(0, s1.basic_vars(row, 0)) = s1.basic_values(row, 0);
		}
		//cout << resolveVarName(s1, s1.basic_vars(row, 0)) << " = " << s1.basic_values(row, 0) << endl;
	}

	// Output all the data
    for(int col = 0; col < int(result.n_cols); ++col)
	{
		cout << resolveVarName(s1, col) << " = " << result(0, col) << endl;
	}

	cout << endl << "Objective Value: " << s1.objvalue << endl << endl;
			
}
//-------------------------------------------------------------------------------------------
// Return the bound that the nonbasic variable is currently resting on, or "active'
//-------------------------------------------------------------------------------------------
double getNonbasicVal(dictionary s1, int index)
{
	if( s1.nonbasic_values(0, index) ) // value is on lower bound
		return s1.nonbasic_upper(0, index);
	else
		return s1.nonbasic_lower(0, index);
}
//-------------------------------------------------------------------------------------------
// Convert a variable name index into a name string
// For example: 0 = x1, 1 = x2, 2 = w1, 3 = w2, 4 = w3
//
// Nonbasic Vars:       0 to (n-1)
// Basic Vars:          n to (n+m-1)
// Auxillary Vars e:    (n+m) to (n+m+n_es-1)
//-------------------------------------------------------------------------------------------
string resolveVarName(dictionary s1, int var_index)
{
    ostringstream strs;
   
	
	// Check if var is nonbasic
	if(var_index < s1.n)
	{
   	    strs << (var_index + 1);
		return "x"+strs.str();		
	}
	// Check if var is basic
	else if( var_index < s1.n + s1.m )
	{
		strs << (var_index - s1.n + 1);		
		return "w"+strs.str();
	}
	// Assume var is auxillary
	else
	{
		strs << (var_index - s1.n - s1.m + 1);		
		return "e"+strs.str();		
	}
}
//-------------------------------------------------------------------------------------------
// Output positive nums less than 10 with extra space
//-------------------------------------------------------------------------------------------
string tabber(double num)
{
    ostringstream strs;

	if(fabs(num) < 0.000000001)
	{
		// assume zero
		return "0  ";
	}
    if(num >= 0 && num < 10)
    {
		strs << num << "  ";
        return strs.str();
    }
    else
    {
		strs << num;		
        return strs.str();
    }
}

//-------------------------------------------------------------------------------------------
// trim from both ends
//-------------------------------------------------------------------------------------------
static inline std::string &trim(std::string &s) {
  return ltrim(rtrim(s));
}

//-------------------------------------------------------------------------------------------
// trim from start
//-------------------------------------------------------------------------------------------
static inline std::string &ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}
//-------------------------------------------------------------------------------------------
// trim from end
//-------------------------------------------------------------------------------------------
static inline std::string &rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}
