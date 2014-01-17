#ifndef HILBERTORDER_HPP
#define HILBERTORDER_HPP 1
#include <vector>
#include <stack>
#include <algorithm>
#include <cmath>
#include "geometry.hpp"
using namespace std;

/*!
\brief Returns the Hilber curve ordering
\param cor The points
\param num The number of points (starting from 0) to calculate the order for
\param innernum The number of points to exclude (starting from 0)
\return The hilbert order indeces
*/
vector<int> HilbertOrder(vector<Vector2D> const& cor,int num,int innernum=0);

#endif // HILBERTORDER_HPP