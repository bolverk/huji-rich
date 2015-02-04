#include "OuterBoundary3D.hpp"
#include "Subcube.hpp"
#include <stdexcept>

using namespace std;

OuterBoundary3D::OuterBoundary3D(OuterBoundary3D::Kinds kind, Vector3D frontUpperRight, Vector3D backLowerLeft) :
	_kind(kind), _frontUpperRight(frontUpperRight), _backLowerLeft(backLowerLeft)
{
	if (kind != PERIODIC && kind != RECTANGULAR)
		throw invalid_argument("OuterBoundary3D kind must be PERIODIC or RECTANGULAR");

	// Make sure frontUpperRight is indeed ahead, above and to the right of _backLowerLeft
	Vector3D diff = _frontUpperRight - _backLowerLeft; // All coordinates must be positive
	if (diff.x <= 0 || diff.y <= 0 || diff.z <= 0)
		throw invalid_argument("In OuterBoundary3D FrontUpperRight must be ahead, above and to the right of BackLowerLeft");
}

OuterBoundary3D::OuterBoundary3D() :
	_kind(RECTANGULAR), _frontUpperRight(Vector3D(1, 1, 1)), _backLowerLeft(Vector3D(0, 0, 0))
{

}

double OuterBoundary3D::distance(const Vector3D &pt, Subcube subcube) const
{
	switch (subcube.NonCenters())
	{
	case 0:
		return 0.0;
	case 1:  // Subcube is one of the 6 cubes sharing a common face
		return distance_face(pt, subcube);
	case 2:  // Subcube is one of the 12 cubes sharing a common edge
		return distance_edge(pt, subcube);
	case 3:  // Subcube is at one of the corners, sharing just a point
		return distance_point(pt, subcube);
	}

	BOOST_ASSERT(false); // How can there be more than 3???
	return -1;
}

//\brief Calculate the distance to a face
//\remark Only one of subcube offset is minus or plus.
//\remark The distance to the face is simply the difference in the axis's coordinate. For example
// for the subcube "+  ", we take _frontUpperRight.x - pt.x
double OuterBoundary3D::distance_face(const Vector3D &pt, Subcube subcube) const
{
	BOOST_ASSERT(subcube.NonCenters() == 1);
	double dist;

	if (subcube[0] == '-')
		dist = pt.x - _backLowerLeft.x;
	else if (subcube[0] == '+')
		dist = pt.x - _frontUpperRight.x;
	else if (subcube[1] == '-')
		dist = pt.y - _backLowerLeft.y;
	else if (subcube[1] == '+')
		dist = pt.y - _frontUpperRight.y;
	else if (subcube[2] == '-')
		dist = pt.z - _backLowerLeft.z;
	else if (subcube[2] == '+')
		dist = pt.z - _frontUpperRight.z;

	return abs(dist);
}

//\brief Calculate the distance between pt to a corner of the cube
//\remark We build the corner from the outer-boundary, and calculate the distance as usual
double OuterBoundary3D::distance_point(const Vector3D &pt, Subcube subcube) const
{
	double cx = subcube[0] == '-' ? _backLowerLeft.x : _frontUpperRight.x;
	double cy = subcube[1] == '-' ? _backLowerLeft.y : _frontUpperRight.y;
	double cz = subcube[2] == '-' ? _backLowerLeft.z : _frontUpperRight.z;

	Vector3D corner(cx, cy, cz);
	return abs(corner - pt);
}

//\brief Calculate the distance between pt to an edge of the cube
//\remark We actually calculate the distance in 2D. For example, for the subcube '++ ', 
// we take the top right edge (X and Y are fixed, Z going along the main subcube). The distance
// between pt and this edge is the distance between (pt.x, pt.y) and (X,Y). This is basically a projection
// to a 2D plane
// \remark To keep the code a bit more concise, we introduce two axis-aliases - U and V, for the two
// non-fixed axes.
double OuterBoundary3D::distance_edge(const Vector3D &pt, Subcube subcube) const
{
	double ptU, ptV, boundaryU, boundaryV;
	if (subcube[0] == Subcube::CENTER)  // X is the same, Y and Z are different
	{
		ptU = pt.y;
		ptV = pt.z;
		boundaryU = subcube[1] == Subcube::MINUS ? _backLowerLeft.y : _frontUpperRight.y;
		boundaryV = subcube[2] == Subcube::MINUS ? _backLowerLeft.z : _frontUpperRight.z;
	}
	else if (subcube[1] == Subcube::CENTER) // Y is the same, X and Z are different
	{
		ptU = pt.x;
		ptV = pt.z;
		boundaryU = subcube[0] == Subcube::MINUS ? _backLowerLeft.x : _frontUpperRight.x;
		boundaryV = subcube[2] == Subcube::MINUS ? _backLowerLeft.z : _frontUpperRight.z;
	}
	else // X and Y are different, Z is the same
	{
		ptU = pt.x;
		ptV = pt.y;
		boundaryU = subcube[0] == Subcube::MINUS ? _backLowerLeft.x : _frontUpperRight.x;
		boundaryV = subcube[1] == Subcube::MINUS ? _backLowerLeft.y : _frontUpperRight.y;
	}

	return sqrt((ptU - boundaryU) * (ptU - boundaryU) + (ptV - boundaryV) * (ptV - boundaryV));
}

