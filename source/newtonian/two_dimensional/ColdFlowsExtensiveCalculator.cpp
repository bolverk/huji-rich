#include "ColdFlowsExtensiveCalculator.hpp"

namespace 
{
	bool NegativeThermalEnergy(Extensive const& cell)
	{
		if (0.50000001*ScalarProd(cell.momentum, cell.momentum) > cell.energy*cell.mass)
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

	bool IsShocked(size_t index, LinearGaussImproved const& interp, ComputationalCell const& cell,Tessellation
		const& tess,EquationOfState const& eos)
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

ColdFlowsExtensiveCalculator::ColdFlowsExtensiveCalculator
(EquationOfState const& eos,
 LinearGaussImproved const& interp): eos_(eos),interp_(interp){}

namespace 
{
	bool bracketed(int low, int arg, int high)
	{
		return arg >= low && high>arg;
	}
}

void ColdFlowsExtensiveCalculator::operator()
(const vector<Extensive>& fluxes,
const PhysicalGeometry& /*pg*/,
const Tessellation& tess,
const double dt,
const CacheData& cd,
const vector<ComputationalCell>& cells,
vector<Extensive>& extensives,
double /*time*/) const
{
	const vector<Edge>& edge_list = tess.getAllEdges();
	for (size_t i = 0; i<edge_list.size(); ++i)
	{
		const Edge& edge = edge_list[i];
		const Extensive delta = dt*cd.areas[i] * fluxes[i];
		if (bracketed(0, edge.neighbors.first, tess.GetPointNo()))
			extensives[static_cast<size_t>(edge.neighbors.first)] -= 	delta;		
		if (bracketed(0, edge.neighbors.second, tess.GetPointNo()))
			extensives[static_cast<size_t>(edge.neighbors.second)] += delta;
	}
	size_t N = static_cast<size_t>(tess.GetPointNo());
	string entropy = "Entropy";
	for (size_t i = 0; i < N; ++i)
	{
		if (!SmallThermalEnergy(extensives[i]))
			continue;
		if (!IsShocked(i, interp_, cells[i],tess,eos_)||NegativeThermalEnergy(extensives[i]))
		{
			const double Ek = ScalarProd(extensives[i].momentum, extensives[i].momentum) / (2 * extensives[i].mass);
			const double density = extensives[i].mass/cd.volumes[i];			
			extensives[i].energy = Ek + eos_.dp2e(density,eos_.sd2p(safe_retrieve(cells[i].tracers, entropy),density))
				*extensives[i].mass;
		}
	}
}
