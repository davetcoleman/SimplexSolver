all: 
	g++ -I armadillo-2.2.4/include -O1 -Wall -g utility.cpp solver.cpp -o solver


# -Wall 	show warnings?
# -01		enable optimization for armadillo
# -g		??
# -o		output to file


# NOT USED:
# -I 		specify the include directory