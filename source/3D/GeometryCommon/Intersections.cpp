#include "Intersections.hpp"


namespace
{
	bool PointInPolygon(Face const& face, Vector3D const& point)
	{
		Vector3D normal = CrossProduct(face.vertices[0] - point, face.vertices[1] - point);
		std::size_t Nloop = face.vertices.size() - 1;
		for (std::size_t i = 0; i < Nloop; ++i)
			if (ScalarProd(CrossProduct(face.vertices[i + 1] - point, face.vertices[(i + 2) % (Nloop + 1)] - point),
				normal) < 0)
				return false;
		return true;
	}

	bool CircleSegmentIntersect(Vector3D const& p0, Vector3D const& p1, Vector3D const& center, double R)
	{
		Vector3D AC = center - p0;
		Vector3D AB = p1 - p0;
		double d = ScalarProd(AC, AB);
		if (d < 0)
		{
			if (abs(AC)>R)
				return false;
			else
				return true;
		}
		double LAB = abs(AB);
		if (d > LAB*LAB)
		{
			if (abs(center - p1) > R)
				return false;
			else
				return true;
		}
		Vector3D closest = p0 + AB*d / (LAB*LAB);
		if (abs(center - closest) > R)
			return false;
		else
			return true;
	}
}

bool FaceSphereIntersections(Face const& face, Sphere const& sphere, Vector3D const& normal)
{	
	// Find plane equation
	double D = -ScalarProd(normal, face.vertices[0]);
	
	// Find intersecting circle
	double d = std::abs(ScalarProd(normal, sphere.center) + D);
	if (d > sphere.radius)
		return false;
	Vector3D circle_center;
	circle_center.x = sphere.center.x - (ScalarProd(normal, sphere.center) + D)*normal.x;
	circle_center.y = sphere.center.y - (ScalarProd(normal, sphere.center) + D)*normal.y;
	circle_center.z = sphere.center.z - (ScalarProd(normal, sphere.center) + D)*normal.z;
	std::size_t Nloop = face.vertices.size();
	double R = sqrt(sphere.radius*sphere.radius - d*d);
	if (PointInPolygon(face, circle_center))
		return true;
	for (std::size_t i = 0; i < Nloop; ++i)
	{
		if (CircleSegmentIntersect(face.vertices[(i + 1) % Nloop], face.vertices[i], circle_center, R))
			return true;
	}
	return false;
}
