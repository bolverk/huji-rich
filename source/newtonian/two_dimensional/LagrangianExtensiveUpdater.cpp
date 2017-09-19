#include "LagrangianExtensiveUpdater.hpp"
#include "../../misc/utils.hpp"

LagrangianExtensiveUpdater::LagrangianExtensiveUpdater(LagrangianFlux const& fc, ExtensiveUpdater const& beu)
	:fc_(fc),beu_(beu) {}

void LagrangianExtensiveUpdater::operator()
(const vector<Extensive>& fluxes,
const PhysicalGeometry& pg,
const Tessellation& tess,
const double dt,
const CacheData& cd,
const vector<ComputationalCell>& cells,
vector<Extensive>& extensives,
double time, TracerStickerNames const& ts) const
{ 
	beu_(fluxes, pg, tess, dt, cd, cells, extensives,time, ts);
	vector<Extensive> extensives_local(extensives);
	Extensive toadd(extensives[0]);
	size_t t_index=static_cast<size_t>(binary_find(ts.tracer_names.begin(),ts.tracer_names.end(),string("Temperature"))-ts.tracer_names.begin());
	assert(t_index<ts.tracer_names.size());
	const vector<Edge>& edge_list = tess.getAllEdges();
	vector<double> dV(tess.GetPointNo(), 0);
	vector<double> dv_ws(dV);
	fc_.edge_vel_.resize(edge_list.size(),0);
	fc_.ws_.resize(edge_list.size(),0);
	for (size_t i = 0; i<edge_list.size(); ++i)
	{
		const Edge& edge = edge_list[i];
		if (edge.neighbors.first < tess.GetPointNo())
		{
			dV[static_cast<size_t>(edge.neighbors.first)] += cd.areas[i] * dt*fc_.edge_vel_[i];
			dv_ws[static_cast<size_t>(edge.neighbors.first)] += cd.areas[i] * dt*fc_.ws_[i];
		}
		if (edge.neighbors.second < tess.GetPointNo())
		{
			dV[static_cast<size_t>(edge.neighbors.second)] -= cd.areas[i] * dt*fc_.edge_vel_[i];
			dv_ws[static_cast<size_t>(edge.neighbors.second)] -= cd.areas[i] * dt*fc_.ws_[i];
		}
	}
	for (size_t i = 0; i < edge_list.size(); ++i)
	{
		double ws = fc_.ws_[i];
		Edge const& edge = tess.GetEdge(static_cast<int>(i));
		if (tess.GetOriginalIndex(edge.neighbors.first) == tess.GetOriginalIndex(edge.neighbors.second))
			continue;
		double L = cd.areas[i];
		if (edge.neighbors.first < tess.GetPointNo())
		{
			ReplaceExtensive(toadd, extensives_local[static_cast<size_t>(edge.neighbors.first)]);
			toadd *= (L*(ws)*dt / (cd.volumes[static_cast<size_t>(edge.neighbors.first)] + 
					       dV[static_cast<size_t>(edge.neighbors.first)] - 
					       dv_ws[static_cast<size_t>(edge.neighbors.first)]));
			extensives[static_cast<size_t>(edge.neighbors.first)] -= toadd; 
		}
		if (edge.neighbors.second < tess.GetPointNo())
		{
			ReplaceExtensive(toadd, extensives_local[static_cast<size_t>(edge.neighbors.second)]);
			toadd *= (L*(ws)*dt / (cd.volumes[static_cast<size_t>(edge.neighbors.second)] + 
					       dV[static_cast<size_t>(edge.neighbors.second)] + 
					       dv_ws[static_cast<size_t>(edge.neighbors.second)]));
			extensives[static_cast<size_t>(edge.neighbors.second)] += toadd;
		}
	}
	fc_.Reset();
}
