#ifndef INTERSECTIONS_HPP
#define INTERSECTIONS_HPP 1

#include "Face.hpp"

struct Sphere
{
	Vector3D center;
	double radius;
};

bool FaceSphereIntersections(Face const& face, Sphere const& sphere);

#endif //INTERSECTIONS_HPP