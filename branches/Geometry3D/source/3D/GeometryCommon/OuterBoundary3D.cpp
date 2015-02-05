#include "OuterBoundary3D.hpp"
#include "Subcube.hpp"
#include <stdexcept>

using namespace std;

OuterBoundary3D::OuterBoundary3D(Vector3D frontUpperRight, Vector3D backLowerLeft) :
	_frontUpperRight(frontUpperRight), _backLowerLeft(backLowerLeft)
{
	// Make sure frontUpperRight is indeed ahead, above and to the right of _backLowerLeft
	Vector3D diff = _frontUpperRight - _backLowerLeft; // All coordinates must be positive
	if (diff.x <= 0 || diff.y <= 0 || diff.z <= 0)
		throw invalid_argument("In OuterBoundary3D FrontUpperRight must be ahead, above and to the right of BackLowerLeft");
}

/*
OuterBoundary3D::OuterBoundary3D() :
	_frontUpperRight(Vector3D(1, 1, 1)), _backLowerLeft(Vector3D(0, 0, 0))
{

} */

double OuterBoundary3D::distance(const Vector3D &pt, Subcube subcube) const
{
	Vector3D vec = vector(pt, subcube);
	return abs(vec);
}

Vector3D OuterBoundary3D::vector(const Vector3D &pt, Subcube subcube) const
{
	switch (subcube.NonCenters())
	{
	case 0:
		return Vector3D(0,0,0);
	case 1:  // Subcube is one of the 6 cubes sharing a common face
		return vector_face(pt, subcube);
	case 2:  // Subcube is one of the 12 cubes sharing a common edge
		return vector_edge(pt, subcube);
	case 3:  // Subcube is at one of the corners, sharing just a point
		return vector_point(pt, subcube);
	}

	BOOST_ASSERT(false); // How can there be more than 3???
	return Vector3D(0,0,0);
}

//\brief Calculate the vector to a face
//\remark Only one of subcube offset is minus or plus.
//\remark The vector to the face runs along ones of the axes. For example, for subcube '-  ' the vector
//is (boundary_x - pt.x, 0, 0) where boundary_x is the X coordinate of the appropriate boundary face
Vector3D OuterBoundary3D::vector_face(const Vector3D &pt, Subcube subcube) const
{
	BOOST_ASSERT(subcube.NonCenters() == 1);

	if (subcube[0] == '-')
		return Vector3D(_backLowerLeft.x - pt.x, 0, 0);
	else if (subcube[0] == '+')
		return Vector3D(_frontUpperRight.x - pt.x, 0, 0);
	else if (subcube[1] == '-')
		return Vector3D(0, _backLowerLeft.y - pt.y, 0);
	else if (subcube[1] == '+')
		return Vector3D(0, _frontUpperRight.y - pt.y, 0);
	else if (subcube[2] == '-')
		return Vector3D(0, 0, _backLowerLeft.z - pt.z);
	else if (subcube[2] == '+')
		return Vector3D(0, 0, _frontUpperRight.z - pt.z);

	BOOST_ASSERT(false);  // We really should be here due to the first assertion
	return Vector3D(0, 0, 0);
}

//\brief Calculate the vector between pt to a corner of the cube
//\remark We build the corner from the outer-boundary
Vector3D OuterBoundary3D::vector_point(const Vector3D &pt, Subcube subcube) const
{
	double cx = subcube[0] == '-' ? _backLowerLeft.x : _frontUpperRight.x;
	double cy = subcube[1] == '-' ? _backLowerLeft.y : _frontUpperRight.y;
	double cz = subcube[2] == '-' ? _backLowerLeft.z : _frontUpperRight.z;

	Vector3D corner(cx, cy, cz);
	return corner - pt;
}

//\brief Calculate the distance between pt to an edge of the cube
//\remark The point on the edge that's closest to pt gets two coordinates from the edge and the third from the point.
// For example, if the subcube is '-+ ', we need to distance from the edge with X=backLowerLeft.x, Y=frontUpperRight.y. Z is pt.Z,
// since this is the closest point on the edge to pt.
Vector3D OuterBoundary3D::vector_edge(const Vector3D &pt, Subcube subcube) const
{
	Vector3D b(pt);
	if (subcube[0] == Subcube::CENTER)  // X is the same, Y and Z are different
	{
		b.y = subcube[1] == Subcube::MINUS ? _backLowerLeft.y : _frontUpperRight.y;
		b.z = subcube[2] == Subcube::MINUS ? _backLowerLeft.z : _frontUpperRight.z;
	}
	else if (subcube[1] == Subcube::CENTER) // Y is the same, X and Z are different
	{
		b.x = subcube[0] == Subcube::MINUS ? _backLowerLeft.x : _frontUpperRight.x;
		b.z = subcube[2] == Subcube::MINUS ? _backLowerLeft.z : _frontUpperRight.z;
	}
	else // X and Y are different, Z is the same
	{
		b.x = subcube[0] == Subcube::MINUS ? _backLowerLeft.x : _frontUpperRight.x;
		b.y = subcube[1] == Subcube::MINUS ? _backLowerLeft.y : _frontUpperRight.y;
	}

	return b - pt;
}


RectangularBoundary3D::RectangularBoundary3D(Vector3D frontUpperRight, Vector3D backLowerLeft)
	: OuterBoundary3D(frontUpperRight, backLowerLeft)
{
}

Vector3D RectangularBoundary3D::ghost(const Vector3D &pt, const Subcube subcube) const
{
	Vector3D toBoundary = vector(pt, subcube);
	return pt + 2 * toBoundary; // Reflect through the boundary
}
