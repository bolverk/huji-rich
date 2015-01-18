#include "OuterBoundary3D.hpp"

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
