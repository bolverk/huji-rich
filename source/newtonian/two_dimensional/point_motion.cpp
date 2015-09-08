#include "point_motion.hpp"

PointMotion::~PointMotion(void) {}

void PointMotion::ApplyFix(Tessellation const& /*tess*/, vector<ComputationalCell> const& /*cells*/, double /*time*/,
	double /*dt*/, vector<Vector2D> & /*velocities*/)
{}
