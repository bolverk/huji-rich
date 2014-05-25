/*! \file square_grid.hpp
  \brief Creates a square cartesian grid
  \author Almog Yalinewich
  \deprecated Replaced by mesh_generator.hpp
 */

#ifndef SQUARE_GRID_HPP
#define SQUARE_GRID_HPP 1

#include <vector>
#include "../../tessellation/geometry.hpp"

using namespace std;

/*! \brief Square cartesian grid
  \param side Side of the square
  \param np Number of points
  \return Square cartesian grid
 */
vector<Vector2D> square_grid(double side, int np);

/*! \brief Shifts all points by a certain offset
  \param original Unperturbed grid
  \param offset Offset vector
  \return Shifted grid
 */
vector<Vector2D> offset_grid(vector<Vector2D> const& original,
			     Vector2D const& offset);

#endif // SQUARE_GRID_HPP
