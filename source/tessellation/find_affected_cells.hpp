/*! \file find_affected_cells.hpp
  \author Almog Yalinewich
  \brief Determines which cells might be affected by boundary conditions
 */

#ifndef FIND_AFFECTED_CELLS_HPP
#define FIND_AFFECTED_CELLS_HPP

#include <vector>
#include "tessellation.hpp"
#include "shape_2d.hpp"

using std::vector;

/*! \brief Recursively finds all cells that intersect a circle
  \param tess Tessellation
  \param index Current cell index
  \param circle Circle
  \param res List of cell indices that intersect the circle
 */
void find_affected_cells(const Tessellation& tess,
				int index,
				const Circle& circle,
				vector<int> & res);

#endif // FIND_AFFECTED_CELLS
