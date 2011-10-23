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
// Read the problem file in .txt format and convert to matrix form
//-------------------------------------------------------------------------------------------
genForm readProblem()
{
	// Create new general form problem
	genForm inputGenForm;

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
    inputGenForm.m = convertCell(cell);

    // Get the number of variables
    getline(lineStream,cell,',');
    inputGenForm.n = convertCell(cell);

    // --------------------------------------------------------------------
    // The second line has the objective coeffients

    std::getline(data,line);
    inputGenForm.c = readRow(inputGenForm.n, line);

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

    inputGenForm.A = A;

    // --------------------------------------------------------------------
    // The fourth to last line is lower boud vector
    std::getline(data,line);
    inputGenForm.a = readRow(inputGenForm.m, line);

    // --------------------------------------------------------------------
    // 3rd to last: uppper bound vector
    std::getline(data,line);
    inputGenForm.b = readRow(inputGenForm.m, line);

    // --------------------------------------------------------------------
    // 2nd to last: var lower bound
    std::getline(data,line);
    inputGenForm.l = readRow(inputGenForm.n, line);

    // --------------------------------------------------------------------
    // last: var upper bounds
    std::getline(data,line);
    inputGenForm.u = readRow(inputGenForm.n, line);

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
void printGenForm(genForm genForm1)
{
    cout << "--------------------------------------" << endl;
    cout << "Number of Contraints: " << genForm1.m << endl;
    cout << "Number of Variables: " << genForm1.n  << endl;
    cout << "Objective Coeffcients ----------------" << endl << genForm1.c;
    cout << "A ------------------------------------" << endl << genForm1.A;
    cout << "Lower Bound Vector -------------------" << endl << genForm1.a;
    cout << "Upper Bound Vector -------------------" << endl << genForm1.b;
    cout << "Lower Bound Variable -----------------" << endl << genForm1.l;
    cout << "Upper Bound Variable -----------------" << endl << genForm1.u;
    cout << "--------------------------------------" << endl << endl;
}
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
void printSimplexStep(simplex s1)
{
    cout << "--------------------------------------------------------------------------" << endl;
    // Print lower bounds
    cout << "l \t \t |   ";
    for(double i = 0; i <= s1.nonbasic_lower.n_rows; ++i)
    {
        // check if the nonbasic variable is currently on this bound
        if(s1.nonbasic_values(0,i) == 0)
	{
	    cout << "\033[1;31m"<< s1.nonbasic_lower(0,i) << "\033[0m \t";
	}
	else
	{
	    cout << s1.nonbasic_lower(0,i) << " \t ";
	}
    }
    cout << endl;

    // Print upper bounds
    cout << "\t u \t |   ";
    for(double i = 0; i <= s1.nonbasic_upper.n_rows; ++i)
    {
        // check if the nonbasic variable is currently on this bound
        if(s1.nonbasic_values(0,i) == 1)
	{
	    cout << "\033[1;31m"<< s1.nonbasic_upper(0,i) << "\033[0m \t";
	}
	else
	{
	    cout << s1.nonbasic_upper(0,i) << " \t ";
	}
    }
    cout << endl;
    cout << "--------------------------------------------------------------------------" << endl;
    
    // Print objective function
    cout << "\t \t |z= ";
    for(double i = 0; i <= s1.nonbasic.n_rows; ++i)
    {
      cout << tabber(s1.nonbasic(0,i)) << " \t";
    }
    cout << endl;
    cout << "--------------------------------------------------------------------------" << endl;

    for(double row = 0; row < s1.basic.n_rows; ++row)
    {
        cout << s1.basic_lower(row,0) << " \t ";
        cout << s1.basic_upper(row,0) << " \t |   ";

	for(double col = 0; col < s1.basic.n_cols; ++col)
	{
  	    cout << tabber(s1.basic(row, col)) << " \t ";
	}

	cout << endl;
    }
    cout << "--------------------------------------------------------------------------" << endl;

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
