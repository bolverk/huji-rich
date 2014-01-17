/*! \file simple_io.hpp
  \brief A collection of simple input / output methods
  \author Almog Yalinewich
*/

#ifndef SIMPLE_IO_HPP
#define SIMPLE_IO_HPP 1

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

/*! \brief Writes a single number to a file
\param num Number
\param fname Name of the file
\param prec Precision (number of digits)
*/
void write_number(double num,
		  string const& fname,
		  int prec=6);

/*! \brief Writes a list of numbers to a file
\param v List of numbers
\param fname Name of the file
\param prec Precision
*/		  
void write_vector(vector<double> const& v,
		  string const& fname,
		  int prec=6);

/*! \brief Writes a list of integers to a file
\param v List of integers
\param fname Name of file
*/
void write_vector(vector<int> const& v,
		  string const& fname);

/*! \brief Reads a single number from a file
\param fname Name of the file
\return The number
*/
double read_number(string const& fname);

/*! \brief Reads a single integer from a file
\param fname Name of the file
\return The integer
*/
int read_int(string const& fname);

#endif // SIMPLE_IO_HPP
