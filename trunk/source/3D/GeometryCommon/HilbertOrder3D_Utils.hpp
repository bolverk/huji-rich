/*! \file HilbertOrder3D_Utils.hpp
\brief Hilbert 3D Order - Utility functions
\author Itai Linial
*/

#ifndef HILBERTORDER3D_UTILS_HPP
#define HILBERTORDER3D_UTILS_HPP 1

#include "Vector3D.hpp"
#include <vector>
#include <algorithm>

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
\param vD_sorted The input points
\param vOut (out) The output points
*/
void FindEqualIndices(vector<size_t> const & vD_sorted, vector<vector<size_t> > & vOut);

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
