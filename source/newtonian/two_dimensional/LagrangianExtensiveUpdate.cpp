#include "LagrangianExtensiveUpdate.hpp"
#include "../../misc/utils.hpp"
namespace
{
	bool bracketed(int low, int arg, int high)
	{
		return arg >= low && high > arg;
	}
}

LagrangianExtensiveUpdate::LagrangianExtensiveUpdate(LagrangianFlux const& lflux, EquationOfState const& eos, GhostPointGenerator
	const& ghost, LinearGaussImproved const& interp) :lflux_(lflux), eos_(eos), ghost_(ghost), interp_(interp) {}

void LagrangianExtensiveUpdate::operator()(const vector<Extensive>& fluxes, const PhysicalGeometry& pg,
	const Tessellation& tess, const double dt, const CacheData& cd, const vector<ComputationalCell>& cells,
	vector<Extensive>& extensives, double time, TracerStickerNames const& tracerstickernames) const
{
	const vector<Edge>& edge_list = tess.getAllEdges();
	Extensive delta = dt*cd.areas[0] * fluxes[0];
	size_t N = tess.GetPointNo();
	std::vector<double> dA(N, 0), dWs(N, 0);
	size_t indexX = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaX")) - tracerstickernames.tracer_names.begin());
	size_t indexY = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaY")) - tracerstickernames.tracer_names.begin());
	assert(indexX < tracerstickernames.tracer_names.size() && indexY < tracerstickernames.tracer_names.size());
	for (size_t i = 0; i < N; ++i)
	{
		extensives[i].tracers[indexX] *= 0.85;
		extensives[i].tracers[indexY] *= 0.85;
	}
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
					&& p_star > 1.2*cells[static_cast<size_t>(edge.neighbors.second)].pressure)
					v_new = std::max(v_new, ScalarProd(cells[static_cast<size_t>(edge.neighbors.second)].velocity,
						normal));
				if (v_new < 0 && tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first)
					&& p_star > 1.2*cells[static_cast<size_t>(edge.neighbors.first)].pressure)
					v_new = std::min(v_new, ScalarProd(cells[static_cast<size_t>(edge.neighbors.first)].velocity,
						normal));
				delta.energy = p_star*v_new*cd.areas[i] * dt;
			}
			else
			{
				if (v_new > 0 && tess.GetOriginalIndex(edge.neighbors.first) != tess.GetOriginalIndex(edge.neighbors.second))
					delta.energy = cd.areas[i] * dt*v_new*cells[static_cast<size_t>(edge.neighbors.first)].pressure;
				else
					if (tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first))
						delta.energy = cd.areas[i] * dt*v_new*cells[static_cast<size_t>(edge.neighbors.second)].pressure;
					else
						delta.energy = 0;
			}
		}
		if (bracketed(0, edge.neighbors.first, tess.GetPointNo()))
		{
			extensives[static_cast<size_t>(edge.neighbors.first)] -= delta;
			extensives[static_cast<size_t>(edge.neighbors.first)].tracers[indexX] -= normal.x*deltaWs;
			extensives[static_cast<size_t>(edge.neighbors.first)].tracers[indexY] -= normal.y*deltaWs;
		}
		if (bracketed(0, edge.neighbors.second, tess.GetPointNo()))
		{
			extensives[static_cast<size_t>(edge.neighbors.second)] += delta;
			extensives[static_cast<size_t>(edge.neighbors.second)].tracers[indexX] -= normal.x*deltaWs;
			extensives[static_cast<size_t>(edge.neighbors.second)].tracers[indexY] -= normal.y*deltaWs;
		}
	}
	// check if cold flows is need
	if (std::binary_search(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(), string("Entropy")))
	{
		ColdFlowsUpdate ceu(eos_, ghost_, interp_);
		for (size_t i = 0; i < N; ++i)
		{
			ceu(fluxes, pg, tess, dt, cd, cells, extensives[i], i, time, tracerstickernames);
		}
	}
}

