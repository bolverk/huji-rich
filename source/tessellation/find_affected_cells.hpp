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

/*! \brief Determines if an edge and a circle intersect
  \param edge Edge
  \param circle Circle
  \return True if both intersect, false otherwise
 */
bool edge_circle_intersect
(const Edge& edge,
 const Circle& circle);

/*! \brief Recursively finds all cells that intersect a circle
  \param tess Tessellation
  \param index Current cell index
  \param circle Circle
  \param res List of cell indices that intersect the circle
 */
void find_affected_cells_recursive(const Tessellation& tess,
				int index,
				const Circle& circle,
				vector<int> & res);

/*! \brief Non recursive version of find affected cells. Only searches one degree of separation
  \param tess Tessellation
  \param index Cell index
  \param circle Circle
  \param vtemp Temperaroy object for not reallocating on heap
  \return Vector with indices of affected cells
 */
vector<int> find_affected_cells
(const Tessellation& tess,
 int index,
 const Circle& circle, vector<int> &vtemp);

#endif // FIND_AFFECTED_CELLS
