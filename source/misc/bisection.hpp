/*! \file bisection.hpp
  \brief Solution of transcendental equations with a single variable using the bisection method
  \author Almog Yalinewich
*/

#ifndef BISECTION_HPP
#define BISECTION_HPP 1

#include "func_1_var.hpp"

/*! \brief Finds an upper limit to an equation
\details Doubles the value until the equation changes sign
\param f One dimensional function
\param xl Lower limit
\return The upper limit
*/
double find_upper_bracket(Func1Var const& f,
			  double xl);

/*! \brief Solves a monotonous transcendental equation using the bisection method
\param f One dimensional function
\param xl Left bound
\param xr Right bound
\param tol Tolerance
\return The solution
*/
double bisection(Func1Var const& f,
		 double xl,
		 double xr,
		 double tol = 1e-6);

#endif // BISECTION_HPP
