//\file Subcube.cpp
//\brief Implementation of the subcube class - representing one of the 27 subcubes
//\author Itay Zandbank

#include "Subcube.h"

using namespace std;

Subcube::Subcube(const char offsets[3])
{
	for (int i = 0; i < 3; i++)
	{
		BOOST_ASSERT(offsets[i] == MINUS || offsets[i] == CENTER || offsets[i] == PLUS);
		_offsets[i] = offsets[i];
	}
}

bool operator<(const Subcube &sc1, const Subcube &sc2)
{
	return sc1.Num() < sc2.Num();
}

bool operator==(const Subcube &sc1, const Subcube &sc2)
{
	return sc1.Num() == sc2.Num();
}

static double distance_face(const Vector3D &pt, const OuterBoundary3D &boundary, Subcube subcube);
static double distance_edge(const Vector3D &pt, const OuterBoundary3D &boundary, Subcube subcube);
static double distance_point(const Vector3D &pt, const OuterBoundary3D &boundary, Subcube subcube);

double distance(const Vector3D &pt, const OuterBoundary3D &boundary, Subcube subcube)
{
	switch (subcube.NonCenters())
	{
	case 0:
		return 0.0;
	case 1:  // Subcube is one of the 6 cubes sharing a common face
		return distance_face(pt, boundary, subcube);
	case 2:  // Subcube is one of the 12 cubes sharing a common edge
		return distance_edge(pt, boundary, subcube);
	case 3:  // Subcube is at one of the corners, sharing just a point
		return distance_point(pt, boundary, subcube);
	}

	BOOST_ASSERT(false); // How can there be more than 3???
	return -1;
}

//\brief Calculate the distance to a face
//\remark Only one of subcube offset is minus or plus.
//\remark The distance to the face is simply the difference in the axis's coordinate. For example
// for the subcube "+  ", we take boundary.FrontUpperRight().x - pt.x
double distance_face(const Vector3D &pt, const OuterBoundary3D &boundary, Subcube subcube)
{
	BOOST_ASSERT(subcube.NonCenters() == 1);
	double dist;

	if (subcube[0] == '-')
		dist = pt.x - boundary.BackLowerLeft().x;
	else if (subcube[0] == '+')
		dist = pt.x - boundary.FrontUpperRight().x;
	else if (subcube[1] == '-')
		dist = pt.y - boundary.BackLowerLeft().y;
	else if (subcube[1] == '+')
		dist = pt.y - boundary.FrontUpperRight().y;
	else if (subcube[2] == '-')
		dist = pt.z - boundary.BackLowerLeft().z;
	else if (subcube[2] == '+')
		dist = pt.z - boundary.FrontUpperRight().z;

	return abs(dist);
}

//\brief Calculate the distance between pt to a corner of the cube
//\remark We build the corner from the outer-boundary, and calculate the distance as usual
double distance_point(const Vector3D &pt, const OuterBoundary3D &boundary, Subcube subcube)
{
	double cx = subcube[0] == '-' ? boundary.BackLowerLeft().x : boundary.FrontUpperRight().x;
	double cy = subcube[1] == '-' ? boundary.BackLowerLeft().y : boundary.FrontUpperRight().y;
	double cz = subcube[2] == '-' ? boundary.BackLowerLeft().z : boundary.FrontUpperRight().z;

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
double distance_edge(const Vector3D &pt, const OuterBoundary3D &boundary, Subcube subcube)
{
	double ptU, ptV, boundaryU, boundaryV;
	if (subcube[0] == Subcube::CENTER)  // X is the same, Y and Z are different
	{
		ptU = pt.y;
		ptV = pt.z;
		boundaryU = subcube[1] == Subcube::MINUS ? boundary.BackLowerLeft().y : boundary.FrontUpperRight().y;
		boundaryV = subcube[2] == Subcube::MINUS ? boundary.BackLowerLeft().z : boundary.FrontUpperRight().z;
	}
	else if (subcube[1] == Subcube::CENTER) // Y is the same, X and Z are different
	{
		ptU = pt.x;
		ptV = pt.z;
		boundaryU = subcube[0] == Subcube::MINUS ? boundary.BackLowerLeft().x : boundary.FrontUpperRight().x;
		boundaryV = subcube[2] == Subcube::MINUS ? boundary.BackLowerLeft().z : boundary.FrontUpperRight().z;
	}
	else // X and Y are different, Z is the same
	{
		ptU = pt.x;
		ptV = pt.y;
		boundaryU = subcube[0] == Subcube::MINUS ? boundary.BackLowerLeft().x : boundary.FrontUpperRight().x;
		boundaryV = subcube[1] == Subcube::MINUS ? boundary.BackLowerLeft().y : boundary.FrontUpperRight().y;
	}

	return sqrt((ptU - boundaryU) * (ptU - boundaryU) + (ptV - boundaryV) * (ptV - boundaryV));
}
