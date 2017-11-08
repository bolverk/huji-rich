/*! \file HilbertProcPositions.hpp
\brief Assign new position to cpu points based on Hilbert curve
\author Elad Steinberg
*/

#include "../3D/GeometryCommon/HilbertOrder3D.hpp"
#include "../3D/GeometryCommon/Tessellation3D.hpp"

#ifndef HILBERTPROC_HPP
#define HILBERTPROC_HPP 1

vector<Vector3D> HilbertProcPositions(Tessellation3D const& tess);

#endif //HILBERTPROC_HPP
