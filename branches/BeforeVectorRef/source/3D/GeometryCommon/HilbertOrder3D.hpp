/*! \file HilbertOrder3D.hpp
\brief Hilbert 3D-space filling curve
\author Itai Linial
*/

// Need t do, add option to define number of points the same for early break, add protection against two points the same

#ifndef HILBERTORDER3D_HPP
#define HILBERTORDER3D_HPP 1
#include <vector>
#include "../../misc/utils.hpp"
#include "Vector3D.hpp"

using namespace std;

/*!
\brief Returns the 3D-Hilbert curve ordering
\param cor The points
\return The Hilbert order indices
*/
vector<size_t> HilbertOrder3D(vector<Vector3D> const& cor);
//void HilbertOrder3D();

#endif // HILBERTORDER_HPP