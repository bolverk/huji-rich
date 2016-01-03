#include "ConditionExtensiveUpdater.hpp"

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
	vector<Extensive>& extensives) const
{
	const vector<Edge>& edge_list = tess.getAllEdges();
	for (size_t i = 0; i<edge_list.size(); ++i)
	{
		const Edge& edge = edge_list[i];
		const Extensive delta = dt*cd.areas[i] * fluxes[i];
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
			if (sequence_[j].first->operator()(i, tess, cells))
			{
				sequence_[j].second->operator()(fluxes, pg, tess, dt, cd, cells, extensives[i],i);
				break;
			}
		}
	}
}

ColdFlowsUpdate::ColdFlowsUpdate(EquationOfState const& eos, LinearGaussImproved const& interp) :
	eos_(eos), interp_(interp) {}

namespace
{
	bool NegativeThermalEnergy(Extensive const& cell)
	{
		if (0.505*ScalarProd(cell.momentum, cell.momentum) > cell.energy*cell.mass)
			return true;
		else
			return false;
	}

	bool SmallThermalEnergy(Extensive const& cell)
	{
		if (0.5*ScalarProd(cell.momentum, cell.momentum) > 0.95*cell.energy*cell.mass)
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
		if (dv < -0.2*eos.dp2c(cell.density, cell.pressure))
			return true;
		const double dp2 = R*R*(interp.GetSlopesUnlimited()[index].xderivative.pressure*
			interp.GetSlopesUnlimited()[index].xderivative.pressure +
			interp.GetSlopesUnlimited()[index].yderivative.pressure*
			interp.GetSlopesUnlimited()[index].yderivative.pressure);
		if (dp2 > 0.1*cell.pressure*cell.pressure)
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
	size_t index)const
{
	if (!SmallThermalEnergy(extensive))
		return;
	string entropy = "Entropy";
	if (!IsShocked(index, interp_, cells[index], tess, eos_) || NegativeThermalEnergy(extensive))
	{
		const double Ek = ScalarProd(extensive.momentum, extensive.momentum) / (2 * extensive.mass);
		const double density = extensive.mass / cd.volumes[index];
		extensive.energy = Ek + eos_.dp2e(density, eos_.sd2p(cells[index].tracers.at(entropy), density))
			*extensive.mass;
	}
}
