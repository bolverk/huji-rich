#include "LagrangianExtensiveUpdater.hpp"
#include "source/misc/utils.hpp"

LagrangianExtensiveUpdater::LagrangianExtensiveUpdater(LagrangianFluxT const& fc):fc_(fc) {}

void LagrangianExtensiveUpdater::operator()
(const vector<Extensive>& fluxes,
const PhysicalGeometry& /*pg*/,
const Tessellation& tess,
const double dt,
const CacheData& cd,
const vector<ComputationalCell>& /*cells*/,
vector<Extensive>& extensives,
double /*time*/, TracerStickerNames const& /*tracerstickernames*/) const
{ 
	const vector<Edge>& edge_list = tess.getAllEdges();
	Extensive delta = dt*cd.areas[0] * fluxes[0];
	vector<double> dV(tess.GetPointNo(), 0);
	fc_.edge_vel_.resize(edge_list.size(),0);
	fc_.ws_.resize(edge_list.size(),0);
	for (size_t i = 0; i<edge_list.size(); ++i)
	{
		const Edge& edge = edge_list[i];
		ReplaceExtensive(delta, fluxes[i]);
		delta *= dt*cd.areas[i];
		if (edge.neighbors.first < tess.GetPointNo())
		{
			extensives[static_cast<size_t>(edge.neighbors.first)] -= delta;
			dV[static_cast<size_t>(edge.neighbors.first)] += cd.areas[i]*dt*fc_.edge_vel_[i]; 
		}
		if (edge.neighbors.second < tess.GetPointNo())
		{
			extensives[static_cast<size_t>(edge.neighbors.second)] += delta;
			dV[static_cast<size_t>(edge.neighbors.second)] -= cd.areas[i]*dt*fc_.edge_vel_[i];
		}
	}
	for (size_t i = 0; i < edge_list.size(); ++i)
	{
		double ws = fc_.ws_[i];
		Edge const& edge = tess.GetEdge(static_cast<int>(i));
		double L = cd.areas[i];
		if (edge.neighbors.first< tess.GetPointNo())
			extensives[static_cast<size_t>(edge.neighbors.first)] -=
				(L*(ws)*dt / (cd.volumes[edge.neighbors.first]+ dV[edge.neighbors.first])) 
				*extensives[static_cast<size_t>(edge.neighbors.first)];
		if (edge.neighbors.second< tess.GetPointNo())
			extensives[static_cast<size_t>(edge.neighbors.second)] +=
				(L*(ws)*dt / (cd.volumes[edge.neighbors.second] + dV[edge.neighbors.second]))
				*extensives[static_cast<size_t>(edge.neighbors.second)];
	}
}
