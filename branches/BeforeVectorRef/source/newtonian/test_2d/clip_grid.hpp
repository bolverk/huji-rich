/*! \file clip_grid.hpp
  \author Almog Yalinewich
  \brief A function that removes all points outside a shape
 */

#ifndef CLIP_GRID_HPP
#define CLIP_GRID_HPP 1

#include <vector>
#include "../../tessellation/geometry.hpp"
#include "../../tessellation/shape_2d.hpp"

using std::vector;

/*! \brief Removes all the points outside a shape
  \param shape A shape
  \param original List of points
  \return List of points inside the shape
 */
vector<Vector2D> clip_grid(const Shape2D& shape,
			   const vector<Vector2D>& original);

#endif
