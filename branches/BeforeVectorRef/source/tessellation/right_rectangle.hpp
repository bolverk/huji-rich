/*! \file right_rectangle.hpp
  \author Almog Yalinewich
  \brief Defines a rectangular shape
 */

#ifndef RECTANGLE_HPP
#define RECTANGLE_HPP 1

#include "shape_2d.hpp"

//! \brief A rectangle aligned with the primary axis (x and y)
class RightRectangle: public Shape2D
{
public:

  /*! \brief Class constructor
    \param lower_left Vertex with lowest values of x and y
    \param upper_right Vertex with highest values of x and y
   */
  RightRectangle(const Vector2D& lower_left,
		 const Vector2D& upper_right);

  bool operator()(const Vector2D& point) const;

private:
  const Vector2D lower_left_;
  const Vector2D upper_right_;
};

#endif // RECTANGLE_HPP
