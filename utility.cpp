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
dictionary readProblem()
{
	// Create new general form problem
	dictionary inputGenForm;

    // Read the actual file
    std::ifstream  data("chapter_example.txt");

    std::string line;

    // --------------------------------------------------------------------
    // The first line should have the number of constraints and variables
    // number constrains:
    std::getline(data,line);
    std::stringstream  lineStream(line);
    std::string        cell;

    // Get the number of constraints
    getline(lineStream,cell,',');
    int constraints_m = convertCell(cell);

    // Get the number of variables
    getline(lineStream,cell,',');
    int variables_n = convertCell(cell);

    // --------------------------------------------------------------------
    // The second line has the objective coeffients

    std::getline(data,line);
    inputGenForm.nonbasic = readRow(variables_n, line);

    // --------------------------------------------------------------------
    // The next m lines have the matrix rows

    mat A(constraints_m,variables_n);    // Set the size of matrix A

    for(int m = 0; m < constraints_m; ++m)
    {
      getline(data,line);
      std::stringstream  lineStream(line);

      for(int n = 0; n < variables_n; ++n)
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
    inputGenForm.basic_lower = trans(readRow(constraints_m, line));

    // --------------------------------------------------------------------
    // 3rd to last: uppper bound vector
    std::getline(data,line);
    inputGenForm.basic_upper = trans(readRow(constraints_m, line));

    // --------------------------------------------------------------------
    // 2nd to last: var lower bound
    std::getline(data,line);
    inputGenForm.nonbasic_lower = readRow(variables_n, line);

    // --------------------------------------------------------------------
    // last: var upper bounds
    std::getline(data,line);
    inputGenForm.nonbasic_upper = readRow(variables_n, line);

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
void printDictionary(dictionary s1, int step)
{
    cout << "Simplex Step " << step << endl;	
    cout << "--------------------------------------------------------------------------" << endl;
	// Print variable locations
	cout << "\t \t \t |   ";
	for(int i = 0; i < int(s1.nonbasic.n_cols); ++i)
	{
		cout << resolveVarName(s1, s1.nonbasic_vars(0, i)) << " \t";
	}
	cout << endl;
	
    // Print lower bounds
    cout << "\tl \t \t |   ";
    for(int i = 0; i <= int(s1.nonbasic_lower.n_rows); ++i)
    {
        // check if the nonbasic variable is currently on this bound
        if(s1.nonbasic_values(0,i) == 0)
	{
	    cout << "\033[1;31m"<< tabber(s1.nonbasic_lower(0,i)) << "\033[0m \t";
	}
	else
	{
	    cout << tabber(s1.nonbasic_lower(0,i)) << "\t";
	}
    }
    cout << endl;

    // Print upper bounds
    cout << "\t\t u \t |   ";
    for(int i = 0; i <= int(s1.nonbasic_upper.n_rows); ++i)
    {
        // check if the nonbasic variable is currently on this bound
        if(s1.nonbasic_values(0,i) == 1)
	{
	    cout << "\033[1;31m"<< tabber(s1.nonbasic_upper(0,i)) << "\033[0m \t";
	}
	else
	{
	    cout << tabber(s1.nonbasic_upper(0,i)) << " \t";
	}
    }
    cout << endl;
    cout << "--------------------------------------------------------------------------" << endl;
    
    // Print objective function
    cout << "\t\t \t |z= ";
    for(int i = 0; i <= int(s1.nonbasic.n_rows); ++i)
    {
      cout << tabber(s1.nonbasic(0,i)) << " \t";
    }
    cout << endl;
    cout << "--------------------------------------------------------------------------" << endl;

    for(double row = 0; row < int(s1.basic.n_rows); ++row)
    {
		cout << resolveVarName(s1, s1.basic_vars(row, 0)) << "\t";
        cout << s1.basic_lower(row,0) << " \t ";
        cout << s1.basic_upper(row,0) << " \t |   ";
	 
		for(int col = 0; col < int(s1.basic.n_cols); ++col)
		{
			cout << tabber(s1.basic(row, col)) << " \t ";
		}

		cout << endl;
    }
    cout << "--------------------------------------------------------------------------" << endl;

}
//-------------------------------------------------------------------------------------------
// Convert a variable name index into a name string
//-------------------------------------------------------------------------------------------
string resolveVarName(dictionary s1, int var_index)
{
    ostringstream strs;
    
	// Check if var is nonbasic or basic
	if(var_index >= int(s1.nonbasic.n_cols))
	{
		// is a basic var
		strs << (var_index - s1.nonbasic.n_cols + 1);		
		return "w"+strs.str();
	}
	else
	{
		// is a non basic var
   	    strs << (var_index + 1);
		return "x"+strs.str();		
	}
}
//-------------------------------------------------------------------------------------------
// Output positive nums less than 10 with extra space
//-------------------------------------------------------------------------------------------
string tabber(double num)
{
    ostringstream strs;
    strs << num;
    
    if(num >= 0 && num < 10)
    {
        return strs.str() + " ";
    }
    else
    {
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
