#ifndef PREDICATES3D_HPP
#define PREDICATES3D_HPP 1

#include <vector>
#include "Vector3D.hpp"
#include <boost/array.hpp>

double orient3d(boost::array<Vector3D,4> const& points);

double insphere(boost::array<Vector3D, 5> const& points);

#endif //PREDICATES3D_HPP