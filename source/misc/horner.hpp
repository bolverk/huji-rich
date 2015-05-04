/*! \file horner.hpp
  \author Almog Yalinewich
  \brief Evaluates a polynomial using Horner's method
 */

#ifndef HORNER_HPP
#define HORNER_HPP 1

#include <vector>

using std::vector;

double horner(const vector<double> coef_list, double x);

#endif // HORNER_HPP
