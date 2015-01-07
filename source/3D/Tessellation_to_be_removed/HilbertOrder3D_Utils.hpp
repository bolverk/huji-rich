/*! \file HilbertOrder3D_Utils.hpp
\brief Hilbert 3D Order - Utility functions
\author Itai Linial
*/

#ifndef HILBERTORDER3D_UTILS_HPP
#define HILBERTORDER3D_UTILS_HPP 1

#include "Vector3D.hpp"
#include <vector>
#include <algorithm>

using namespace std;

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
\brief Scale a vector of 3D points to the unit-cube
\param vPointsIn The input points
\param vPointsOut (out) The output points
*/
void FindEqualIndices(vector<unsigned long long int> const & vD_sorted, vector<vector<size_t> > & vOut);

// Return indices order after sorting of the values vector:
template <typename T>
vector<size_t> ordered(vector<T> const& values) {
	vector<size_t> indices(values.size());

	for (size_t ii = 0; ii < indices.size(); ++ii)
	{
		indices[ii] = ii;
	}

	sort(
	     //		begin(indices), end(indices),
	     indices.begin(),indices.end(),
		[&](size_t a, size_t b) { return values[a] < values[b]; }
	);
	return indices;
}

// Reorder a vector according to an index vector (obtained from the 'ordered' function)
template< class T >
void reorder(vector<T> & v, vector<size_t> const & order)
{
	vector<T> vCopy = v;
	for (size_t ii = 0; ii < v.size(); ++ii)
	{
		v[ii] = vCopy[order[ii]];
	}
	return;
}

#endif // HILBERTORDER3D_UTILS_HPP
