#include "eulerian_3d.hpp"

Eulerian3D::Eulerian3D(void) {}

void Eulerian3D::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& /*cells*/,
	double /*time*/, TracerStickerNames const& /*tracerstickernames*/, vector<Vector3D> &res) const
{
	res.resize(tess.GetPointNo(), Vector3D());
	return;
}
