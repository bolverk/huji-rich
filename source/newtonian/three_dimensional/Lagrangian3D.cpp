#include "Lagrangian3D.hpp"

void Lagrangian3D::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
	double /*time*/, TracerStickerNames const& /*tracerstickernames*/, vector<Vector3D> &res) const
{
	res.resize(tess.GetPointNo());
	size_t N = res.size();
	for (size_t i = 0; i<N; ++i)
		res[i] = cells[i].velocity;
}
