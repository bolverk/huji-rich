/*! \file triangle_area.hpp
  \author Almog Yalinewich
  \brief Calculates the area of a triangle
*/

#ifndef TRIANGLE_AREA_HPP
#define TRIANGLE_AREA_HPP 1

#include "geometry.hpp"

/*! \brief Calculates the area of a triangle
  \param p1 Vertex #1
  \param p2 Vertex #2
  \param p3 Vertex #3
  \return Area of the triangle
 */
double calc_triangle_area(const Vector2D& p1,
			  const Vector2D& p2,
			  const Vector2D& p3);

#endif // TRIANGLE_AREA_HPP
