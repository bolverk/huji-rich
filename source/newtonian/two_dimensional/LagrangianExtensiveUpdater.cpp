#include "LagrangianExtensiveUpdater.hpp"
#include "../../misc/utils.hpp"

LagrangianExtensiveUpdater::LagrangianExtensiveUpdater(LagrangianFlux const& fc, ExtensiveUpdater const& beu,
	EquationOfState const& eos)
	:fc_(fc),beu_(beu),eos_(eos) {}

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
	size_t Npoints = static_cast<size_t>(tess.GetPointNo());
	vector<double> R(Npoints, 0);
	Extensive toadd(extensives[0]);
	size_t t_index=static_cast<size_t>(binary_find(ts.tracer_names.begin(),ts.tracer_names.end(),string("Temperature"))-ts.tracer_names.begin());
	assert(t_index<ts.tracer_names.size());
	for (size_t i = 0; i < Npoints; ++i)
		R[i] = tess.GetWidth(static_cast<int>(i));
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
		double Cs = 0;
		Vector2D n = normalize(tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first));
		Cs = eos_.dp2c(cells[static_cast<size_t>(edge.neighbors.first)].density,
			cells[static_cast<size_t>(edge.neighbors.first)].pressure,cells[static_cast<size_t>(edge.neighbors.first)].tracers,ts.tracer_names) 
			+ std::abs(fc_.edge_vel_[i] - ScalarProd(n,cells[static_cast<size_t>(edge.neighbors.first)].velocity));
		Cs = std::max(Cs,eos_.dp2c(cells[static_cast<size_t>(edge.neighbors.second)].density,
			cells[static_cast<size_t>(edge.neighbors.second)].pressure,cells[static_cast<size_t>(edge.neighbors.first)].tracers,ts.tracer_names) 
		+ std::abs(fc_.edge_vel_[i]	- ScalarProd(n, cells[static_cast<size_t>(edge.neighbors.second)].velocity)));
		double density_ratio = cells[static_cast<size_t>(edge.neighbors.first)].density /
			cells[static_cast<size_t>(edge.neighbors.second)].density;
		density_ratio = std::max(density_ratio, 1.0 / density_ratio);
		double p_ratio = cells[static_cast<size_t>(edge.neighbors.first)].pressure /
			cells[static_cast<size_t>(edge.neighbors.second)].pressure;
		p_ratio = std::max(p_ratio, 1.0 / p_ratio);
		double L = cd.areas[i];
		if (edge.neighbors.first < tess.GetPointNo() && cells[static_cast<size_t>(edge.neighbors.first)].stickers[t_index]<2.5e8)
		{
			ReplaceExtensive(toadd, extensives[static_cast<size_t>(edge.neighbors.first)]);
			toadd *= (L*(ws)*dt / (cd.volumes[edge.neighbors.first] + dV[edge.neighbors.first] - dv_ws[edge.neighbors.first]));
			extensives[static_cast<size_t>(edge.neighbors.first)] -= toadd; 
		}
		if (edge.neighbors.second < tess.GetPointNo() && cells[static_cast<size_t>(edge.neighbors.second)].stickers[t_index]<2.5e8)
		{
			ReplaceExtensive(toadd, extensives[static_cast<size_t>(edge.neighbors.second)]);
			toadd *= (L*(ws)*dt / (cd.volumes[edge.neighbors.second] + dV[edge.neighbors.second] + dv_ws[edge.neighbors.second]));
			extensives[static_cast<size_t>(edge.neighbors.second)] += toadd;
		}
	}
	fc_.Reset();
}
