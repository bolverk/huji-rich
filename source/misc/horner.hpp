/*! \file horner.hpp
  \author Almog Yalinewich
  \brief Evaluates a polynomial using Horner's method
 */

#ifndef HORNER_HPP
#define HORNER_HPP 1

#include <vector>

using std::vector;

/*! \brief Evaluates a polynomial using the horner method
  \param coef_list List of coefficients
  \param x Abscissa
  \return Ordinate
 */
double horner(const vector<double> coef_list, double x);

#endif // HORNER_HPP
