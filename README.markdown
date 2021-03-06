LP Solver for General Simplex Problems
==========
* Author:  Dave Coleman <davetcoleman@gmail.com>
* License: GNU General Public License, version 3 (GPL-3.0)
* Date: 10/20/2011

A very simple implementation for a class project in C++


INSTALL NOTES
---------

Tested on Ubuntu 11.04 and Mac OSX
       
       cd into directory
       make
       ./solver FILENAME


or, to run the 5 benchmarking tests in the /tests folder just run:
    
	./solver

If this fails to run you probably need the Armadillo library:

   	sudo apt-get install libarmadillo0


IMPLEMENTATION
---------

I have implemented a general form simplex solver for linear programs. Problems are initialized using an intilization phase as described in Sankaranarayanan's General Initialization Phase notes. The primal simplex method for general LP problems is implmented as described in Chapter 9 of Vanderbei's Linear Programming book.

The program is written in C++ but uses the Armadillo matricies libraries for more dynamic and optimized matrix operations.

To ensure termination I used the least subscript/Bland's rule after the program hits the 50th step. The program aborts if it runs for over 1000 pivots as a debugging feature.

5 examples where testing on this solver:
 - Example from Chapter 9 of Vanderbei's Linear Programming Book
 - Example from Sankaranarayanan's General Intialization PDF notes
 - Example from Sankaranarayanan's Matlab code, #1
 - Example from Sankaranarayanan's Lecture Notes 8, Slide 25
 - Example from Chvatal Chapter 3, page 39

The objective functions were solved using glpk solver for comparison.

For benchmarking each of the 5 tests were run 100,000 times and the total time was added up. The program was then optimized to lower the average running time.

BENCHMARKING 500,000 SOLUTIONS

2.97919	     Solver's initial run time after it was finished
0.124841     Changed matricies to reference variables

95.81% speed up, not bad. I'm tired...


