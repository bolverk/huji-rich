#include "LagrangianExtensiveUpdate.hpp"
#include "../../misc/utils.hpp"
namespace
{
	bool bracketed(int low, int arg, int high)
	{
		return arg >= low && high > arg;
	}
}

LagrangianExtensiveUpdate::LagrangianExtensiveUpdate(LagrangianFlux const& lflux, GhostPointGenerator
	const& ghost) :lflux_(lflux), ghost_(ghost) {}

void LagrangianExtensiveUpdate::operator()(const vector<Extensive>& fluxes, const PhysicalGeometry& /*pg*/,
	const Tessellation& tess, const double dt, const CacheData& cd, const vector<ComputationalCell>& cells,
	vector<Extensive>& extensives, double time, TracerStickerNames const& tracerstickernames) const
{
	const vector<Edge>& edge_list = tess.getAllEdges();
	Extensive delta = dt*cd.areas[0] * fluxes[0];
	size_t N = static_cast<size_t>(tess.GetPointNo());
	std::vector<double> dA(N, 0), dWs(N, 0);
	size_t indexX = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaX")) - tracerstickernames.tracer_names.begin());
	size_t indexY = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaY")) - tracerstickernames.tracer_names.begin());
	//assert(indexX < tracerstickernames.tracer_names.size() && indexY < tracerstickernames.tracer_names.size());
	bool area_calc = false;
	if (indexX < tracerstickernames.tracer_names.size())
	{
		area_calc = true;
		for (size_t i = 0; i < N; ++i)
		{
			extensives[i].tracers[indexX] *= 0.85;
			extensives[i].tracers[indexY] *= 0.85;
		}
	}

	boost::container::flat_map<size_t,ComputationalCell> ghosts = ghost_(tess, cells, time, tracerstickernames);
	vector<ComputationalCell> newcells(cells);
	newcells.resize(static_cast<size_t>(tess.GetTotalPointNumber()));
	for (boost::container::flat_map<size_t, ComputationalCell>::iterator it = ghosts.begin(); it != ghosts.end(); ++it)
		newcells[it->first] = it->second;

	for (size_t i = 0; i < edge_list.size(); ++i)
	{
		const Edge& edge = edge_list[i];
		ReplaceExtensive(delta, fluxes[i]);
		delta *= dt*cd.areas[i];
		double deltaWs = -lflux_.ws_[i] * cd.areas[i] * dt;
		Vector2D normal = normalize(tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first));
		if (lflux_.Lag_calc_[i])
		{
			double p_star = ScalarProd(fluxes[i].momentum, normal);
			double v_star = fluxes[i].energy / p_star;
			double v_new = (v_star - lflux_.ws_[i]);
			if (v_new*v_star > 0)
			{
				if (v_new > 0 && tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first)
					&& p_star > 1.2*newcells[static_cast<size_t>(edge.neighbors.second)].pressure)
					v_new = std::max(v_new, ScalarProd(newcells[static_cast<size_t>(edge.neighbors.second)].velocity,
						normal));
				if (v_new < 0 && tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first)
					&& p_star > 1.2*newcells[static_cast<size_t>(edge.neighbors.first)].pressure)
					v_new = std::min(v_new, ScalarProd(newcells[static_cast<size_t>(edge.neighbors.first)].velocity,
						normal));
				delta.energy = p_star*v_new*cd.areas[i] * dt;
			}
			else
			{
				if (v_new > 0 && tess.GetOriginalIndex(edge.neighbors.first) != tess.GetOriginalIndex(edge.neighbors.second))
					delta.energy = cd.areas[i] * dt*v_new*newcells[static_cast<size_t>(edge.neighbors.first)].pressure;
				else
					if (tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first))
						delta.energy = cd.areas[i] * dt*v_new*newcells[static_cast<size_t>(edge.neighbors.second)].pressure;
					else
						delta.energy = 0;
			}
		}
		if (bracketed(0, edge.neighbors.first, tess.GetPointNo()))
		{
			extensives[static_cast<size_t>(edge.neighbors.first)] -= delta;
			if (area_calc)
			{
				extensives[static_cast<size_t>(edge.neighbors.first)].tracers[indexX] -= normal.x*deltaWs;
				extensives[static_cast<size_t>(edge.neighbors.first)].tracers[indexY] -= normal.y*deltaWs;
			}
		}
		if (bracketed(0, edge.neighbors.second, tess.GetPointNo()))
		{
			extensives[static_cast<size_t>(edge.neighbors.second)] += delta;
			if (area_calc)
			{
				extensives[static_cast<size_t>(edge.neighbors.second)].tracers[indexX] -= normal.x*deltaWs;
				extensives[static_cast<size_t>(edge.neighbors.second)].tracers[indexY] -= normal.y*deltaWs;
			}
		}
	}
}

