/*! \file HilbertOrder.hpp
  \brief Hilbert space filling curve
  \author Elad Steinberg
 */

#ifndef HILBERTORDER_HPP
#define HILBERTORDER_HPP 1
#include <vector>
#include <stack>
#include <algorithm>
#include <cmath>
#include "geometry.hpp"
#include "../misc/universal_error.hpp"

using std::vector;
using std::pair;
using std::stack;
using std::max;
using std::min;

/*!
\brief Returns the Hilber curve ordering
\param cor The points
\param num The number of points (starting from 0) to calculate the order for
\param innernum The number of points to exclude (starting from 0)
\return The hilbert order indeces
*/
vector<int> HilbertOrder(vector<Vector2D> const& cor,int num,int innernum=0);

#endif // HILBERTORDER_HPP
