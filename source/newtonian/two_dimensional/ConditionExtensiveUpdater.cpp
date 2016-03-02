#include "ConditionExtensiveUpdater.hpp"
#include "../../misc/utils.hpp"

namespace
{
	bool bracketed(int low, int arg, int high)
	{
		return arg >= low && high>arg;
	}
}

ConditionExtensiveUpdater::Condition::~Condition() {}

ConditionExtensiveUpdater::Action::~Action() {}

ConditionExtensiveUpdater::~ConditionExtensiveUpdater() {}

ConditionExtensiveUpdater::ConditionExtensiveUpdater(const vector<pair<const Condition*, const Action*> >& sequence) :
	sequence_(sequence) {}

void ConditionExtensiveUpdater::operator()(const vector<Extensive>& fluxes,
	const PhysicalGeometry& pg,
	const Tessellation& tess,
	const double dt,
	const CacheData& cd,
	const vector<ComputationalCell>& cells,
	vector<Extensive>& extensives,
	double time,
	TracerStickerNames const& tracerstickernames) const
{
	const vector<Edge>& edge_list = tess.getAllEdges();
	Extensive delta(extensives[0]);
	for (size_t i = 0; i<edge_list.size(); ++i)
	{
		const Edge& edge = edge_list[i];
		ReplaceExtensive(delta, fluxes[i]);
		delta *= dt*cd.areas[i];
		if (bracketed(0, edge.neighbors.first, tess.GetPointNo()))
			extensives[static_cast<size_t>(edge.neighbors.first)] -= delta;
		if (bracketed(0, edge.neighbors.second, tess.GetPointNo()))
			extensives[static_cast<size_t>(edge.neighbors.second)] += delta;
	}
	size_t n = static_cast<size_t>(tess.GetPointNo());
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < sequence_.size(); ++j)
		{
			if (sequence_[j].first->operator()(i, tess, cells,time,tracerstickernames))
			{
				sequence_[j].second->operator()(fluxes, pg, tess, dt, cd, cells, extensives[i],i,time,tracerstickernames);
				break;
			}
		}
	}
}

ColdFlowsUpdate::ColdFlowsUpdate(EquationOfState const& eos, LinearGaussImproved const& interp) :
	eos_(eos), interp_(interp),entropy_index_(-1) {}

namespace
{
	bool NegativeThermalEnergy(Extensive const& cell)
	{
		if (0.5000001*ScalarProd(cell.momentum, cell.momentum) > cell.energy*cell.mass)
			return true;
		else
			return false;
	}

	bool IsShocked(size_t index, LinearGaussImproved const& interp, ComputationalCell const& cell, Tessellation
		const& tess, EquationOfState const& eos)
	{
		const double R = tess.GetWidth(static_cast<int>(index));
		const double dv = R*(interp.GetSlopesUnlimited()[index].xderivative.velocity.x +
			interp.GetSlopesUnlimited()[index].yderivative.velocity.y);
		if (dv < -0.5*eos.dp2c(cell.density, cell.pressure))
			return true;
		const double dp2 = R*R*(interp.GetSlopesUnlimited()[index].xderivative.pressure*
			interp.GetSlopesUnlimited()[index].xderivative.pressure +
			interp.GetSlopesUnlimited()[index].yderivative.pressure*
			interp.GetSlopesUnlimited()[index].yderivative.pressure);
		if (dp2 > 0.2*cell.pressure*cell.pressure)
			return true;
		return false;
	}
}

void ColdFlowsUpdate::operator()
(const vector<Extensive>& /*fluxes*/,
	const PhysicalGeometry& /*pg*/,
	const Tessellation& tess,
	const double /*dt*/,
	const CacheData& cd,
	const vector<ComputationalCell>& cells,
	Extensive& extensive,
	size_t index,
	double /*time*/,
	TracerStickerNames const& ts)const
{
	if (!NegativeThermalEnergy(extensive))
		return;
	if(entropy_index_<0)
		entropy_index_ = static_cast<int>(lower_bound(ts.tracer_names.begin(),ts.tracer_names.end(),string("Entropy")) - ts.tracer_names.begin());	
	size_t e_index=static_cast<size_t>(entropy_index_);
	assert(e_index<extensive.tracers.size());
	if (!IsShocked(index, interp_, cells[index], tess, eos_) || NegativeThermalEnergy(extensive))
	{
		const double Ek = ScalarProd(extensive.momentum, extensive.momentum) / (2 * extensive.mass);
		const double density = extensive.mass / cd.volumes[index];
		extensive.energy = Ek + eos_.dp2e(density, eos_.sd2p(cells[index].tracers[e_index], density))
			*extensive.mass;
	}
}
