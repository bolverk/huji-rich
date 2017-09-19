/*! \file HilbertOrder3D_Utils.hpp
\brief Hilbert 3D Order - Utility functions
\author Itai Linial
*/

#ifndef HILBERTORDER3D_UTILS_HPP
#define HILBERTORDER3D_UTILS_HPP 1

#include "Vector3D.hpp"
#include <vector>
#include <algorithm>
#include "../../misc/utils.hpp"

using std::vector;

/*!
\brief Estimate the number of iterations required in the Hilbert Curve, according to the number of points
\param cor The points
\return The estimated number of required iterations
*/
int EstimateHilbertIterationNum(vector<Vector3D> const& cor);

/*!
\brief Scale a vector of 3D points to the unit-cube
\param vPointsIn The input points
\param vPointsOut (out) The output points
*/
void AdjustPoints(vector<Vector3D> const & vPointsIn, vector<Vector3D> & vPointsOut);

/*!
\brief Find points with same indeces
\param vD_sorted The input points, sorted
\param vOut (out) The list of points with same indeces
*/
void FindEqualIndices(vector<unsigned long long int> const & vD_sorted, vector<vector<std::size_t> > & vOut);

// Return indices order after sorting of the values vector:
template <typename T>
void ordered(vector<T> const& values, vector<std::size_t> &indices)
{
	indices.resize(values.size());
	sort_index(values,indices);
}

// Reorder a vector according to an index vector (obtained from the 'ordered' function)
template< class T >
void reorder(vector<T> & v, vector<std::size_t> const & order)
{
	vector<T> vCopy = v;
	for (std::size_t ii = 0; ii < v.size(); ++ii)
	{
		v[ii] = vCopy[order[ii]];
	}
	return;
}

#endif // HILBERTORDER3D_UTILS_HPP