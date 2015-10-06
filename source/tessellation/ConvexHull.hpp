/*! \file ConvexHull.hpp
  \brief Calculates the convex hull from a set of points
  \author Elad Steinberg
 */

#ifndef CONVEXHULL_HPP
#define CONVEXHULL_HPP 1

#include "../misc/utils.hpp"
#include "tessellation.hpp"
#include <cmath>

/*! \brief Returns the ConvexHull for a set of points
  \param result The set of convex hull points
  \param tess The tessellation
  \param index The index of the cell for which to calculate the convex hull
*/
void ConvexHull(vector<Vector2D> &result,Tessellation const& tess,int index);
/*! \brief Returns the ConvexHull of the edges of a cell
  \param result The indeces of convex hull edges
  \param tess The tessellation
  \param index The index of the cell for which to calculate the convex hull
*/
void ConvexEdges(vector<int> &result,Tessellation const& tess,int index);
#endif //CONVEXHULL_HPP
