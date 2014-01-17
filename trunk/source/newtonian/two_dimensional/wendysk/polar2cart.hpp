/*! \file polar2cart.hpp
  \brief Converts polar coordinates to cartesian coordinates
  \author Almog Yalinewich
 */

#ifndef POLAR2CART_HPP
#define POLAR2CART_HPP 1

#include <cmath>
#include "../../../tessellation/geometry.hpp"

/*! \brief Converts polar coordinates to a cartesian vector
  \param radius Radius
  \param angle Angle (relative to the x axis)
  \return Two dimensional cartesian vector
 */
Vector2D polar2cart(double radius,
		    double angle);

#endif // POLAR2CART_HPP

