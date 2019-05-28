#include "zero_force.hpp"

vector<Extensive> ZeroForce::operator()
(const Tessellation& tess,
	const PhysicalGeometry& /*pg*/,
	const CacheData& /*cd*/,
 const vector<ComputationalCell>& /*cells*/,
	const vector<Extensive>& /*fluxes*/,
	const vector<Vector2D>& /*point_velocities*/,
	const double /*t*/,
	TracerStickerNames const& /*tracerstickersnames*/) const
{
	vector<Extensive> res(static_cast<size_t>(tess.GetPointNo()));
//	size_t N = res.size();
	//for (size_t i = 0; i < N; ++i)
		//res[i].tracers.resize(cells[0].tracers.size(), 0);
	return res;
}
