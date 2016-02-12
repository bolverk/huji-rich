#include "point_motion.hpp"

PointMotion::~PointMotion(void) {}

vector<Vector2D> PointMotion::ApplyFix(Tessellation const& /*tess*/, vector<ComputationalCell> const& /*cells*/, double /*time*/,
	double /*dt*/, vector<Vector2D> const& velocities, TracerStickerNames const&
	/*tracerstickernames*/)const
{
	return velocities;
}
