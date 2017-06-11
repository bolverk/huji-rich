#ifndef INTERSECTIONS_HPP
#define INTERSECTIONS_HPP 1

#include "Face.hpp"

struct Sphere
{
	Vector3D center;
	double radius;
	Sphere():center(Vector3D()),radius(0){}
};

bool FaceSphereIntersections(Face const& face, Sphere const& sphere,Vector3D const& normal);

#endif //INTERSECTIONS_HPP