#include "point_motion_3d.hpp"

PointMotion3D::~PointMotion3D(void) {}

void PointMotion3D::ApplyFix(Tessellation3D const& /*tess*/, vector<ComputationalCell3D> const& /*cells*/, double /*time*/,
	double /*dt*/, vector<Vector3D> &/*velocities*/, TracerStickerNames const& /*tracerstickernames*/)const
{
	return;
}

