#include "ConditionExtensiveUpdater.hpp"
#include "../../misc/utils.hpp"
#include <iostream>
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
			if (sequence_[j].first->operator()(i, tess, cells, time, tracerstickernames))
			{
				sequence_[j].second->operator()(fluxes, pg, tess, dt, cd, cells, extensives[i], i, time, tracerstickernames);
				break;
			}
		}
	}
}

ColdFlowsUpdate::ColdFlowsUpdate(EquationOfState const& eos, GhostPointGenerator const& ghost,LinearGaussImproved
	const& interp) :
	eos_(eos), ghost_(ghost),interp_(interp),lasttime_(0), dt_(0), entropy_index_(-1), ghost_cells_(boost::container::
		flat_map<size_t, ComputationalCell>()) {}

namespace
{
	Vector2D GetTemperatureGrad(size_t index,LinearGaussImproved const& interp,ComputationalCell const& cell)
	{
		Vector2D res;
		Slope slope=interp.GetSlopesUnlimited()[index];
		res.x = slope.xderivative.pressure / cell.density - slope.xderivative.density*cell.pressure / (cell.density*cell.density);
		res.y = slope.yderivative.pressure / cell.density - slope.yderivative.density*cell.pressure / (cell.density*cell.density);
		return res;
	}

	bool NegativeVelocityDivergence(size_t index, LinearGaussImproved const& interp,double cs,double R)
	{
		Slope slope = interp.GetSlopesUnlimited()[index];
		return (slope.xderivative.velocity.x + slope.yderivative.velocity.y)*R < -0.2*cs;
	}

/*	bool PositiveTemperatureDensityGrad(size_t index, LinearGaussImproved const& interp,Vector2D const& Tgrad)
	{
		Slope slope = interp.GetSlopesUnlimited()[index];
		return (Tgrad.x*slope.xderivative.density + Tgrad.y*slope.yderivative.density)>0;
	}
*/
	bool BigJump(Vector2D const& Tgrad, size_t index,Tessellation const& tess,vector<ComputationalCell> const& cells,
		boost::container::flat_map<size_t, ComputationalCell> ghost_cells)
	{
		vector<int> neigh = tess.GetNeighbors(static_cast<int>(index));
		size_t n = neigh.size();
		int N = tess.GetPointNo();
		double maxT=0, maxP=0;
		const Vector2D & point = tess.GetMeshPoint(static_cast<int>(index));
		for (size_t i = 0; i < n; ++i)
		{
			if (ScalarProd(point - tess.GetMeshPoint(neigh[i]), Tgrad)<0)
			{
				if (neigh[i] < N || (tess.GetOriginalIndex(neigh[i])!=static_cast<int>(index))) 
				{
					maxT = std::max(maxT, cells[static_cast<size_t>(neigh[i])].pressure /
						cells[static_cast<size_t>(neigh[i])].density);
					maxP = std::max(maxP, cells[static_cast<size_t>(neigh[i])].pressure);
				}
				else
				{
					const boost::container::flat_map<size_t, ComputationalCell>::const_iterator it =
						ghost_cells.find(static_cast<size_t>(neigh[i]));
					assert(it != ghost_cells.end());
					maxT = std::max(maxT, it->second.pressure / it->second.density);
					maxP = std::max(maxP, it->second.pressure);
				}
			}
		}
		return (log(maxP) - log(cells[index].pressure))>0.2 || 
			(log(maxT)-log(cells[index].pressure/cells[index].density))>0.1;
	}

	bool NegativeThermalEnergy(Extensive const& cell)
	{
		if (0.50000001*ScalarProd(cell.momentum, cell.momentum) > cell.energy*cell.mass)
			return true;
		else
			return false;
	}

	bool SmallThermalEnergy(Extensive const& cell)
	{
		if (0.51*ScalarProd(cell.momentum, cell.momentum) > cell.energy*cell.mass)
			return true;
		else
			return false;
	}

	void EntropyFix(Extensive &extensive, double vol, EquationOfState const& eos, size_t entropy_index)
	{
		const double Ek = ScalarProd(extensive.momentum, extensive.momentum) / (2 * extensive.mass);
		const double density = extensive.mass / vol;
		extensive.energy = Ek + eos.dp2e(density, eos.sd2p(extensive.tracers[entropy_index] / extensive.mass, density))
			*extensive.mass;
	}

	/*bool EntropyIncrease(ComputationalCell const& oldcell, size_t entropyindex, double new_entropy)
	{
		return oldcell.tracers[entropyindex] * 1.00067 < new_entropy;
	}
*/
	double NewPressure(Extensive const& extensive, EquationOfState const& eos, double new_d)
	{
		double new_p = extensive.energy - ScalarProd(extensive.momentum, extensive.momentum)*0.5 / extensive.mass;
		new_p /= extensive.mass;
		if (new_p < 0)
			return -1;
		return eos.de2p(new_d, new_p);
	}

/*	bool IsShocked(size_t index, ComputationalCell const& cell, Tessellation
		const& tess, vector<int> &temp, vector<ComputationalCell> const& cells,
		boost::container::flat_map<size_t, ComputationalCell> ghost_cells, double new_p, double new_d)
	{
		if (new_p*cell.density<1.1*new_d*cell.pressure)
			return false;
		temp = tess.GetCellEdges(static_cast<int>(index));
		size_t N = temp.size();
		Vector2D normal;
		size_t other_cell = 0;
		for (size_t i = 0; i < N; ++i)
		{
			double other_p = 0;
			Edge const& edge = tess.GetEdge(temp[i]);
			if (edge.neighbors.first == static_cast<int>(index))
			{
				normal = normalize(tess.GetMeshPoint(edge.neighbors.first) - tess.GetMeshPoint(edge.neighbors.second));
				other_cell = static_cast<size_t>(edge.neighbors.second);
			}
			else
			{
				normal = normalize(tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first));
				other_cell = static_cast<size_t>(edge.neighbors.first);
			}
			const boost::container::flat_map<size_t, ComputationalCell>::const_iterator it =
				ghost_cells.find(other_cell);
			const ComputationalCell* cell_other;
			if (it == ghost_cells.end())
				cell_other = &cells.at(other_cell);
			else
				cell_other = &it->second;
			double inout = ScalarProd(normal, cell_other->velocity - cell.velocity) > 0 ? 1 : -1;
			other_p = cell_other->pressure + (inout>0 ? 3 : 1)*cell_other->density*ScalarProd(cell.velocity - cell_other->velocity,
				normal)*ScalarProd(cell.velocity - cell_other->velocity, normal)*inout;
			if (other_p > 1.1*cell.pressure)
				return true;
		}
		return false;
	}*/
}

void ColdFlowsUpdate::operator()
(const vector<Extensive>& /*fluxes*/,
	const PhysicalGeometry& /*pg*/,
	const Tessellation& tess,
	const double dt,
	const CacheData& cd,
	const vector<ComputationalCell>& cells,
	Extensive& extensive,
	size_t index,
	double time,
	TracerStickerNames const& ts)const
{
	if (entropy_index_ < 0)
	{
		entropy_index_ = static_cast<int>(lower_bound(ts.tracer_names.begin(), ts.tracer_names.end(), string("Entropy")) - ts.tracer_names.begin());
		lasttime_ = time;
		ghost_cells_ = ghost_.operator()(tess, cells, time, ts);
		dt_ = dt;
	}
	if (lasttime_ < time || dt < dt_ || dt > dt_)
	{
		lasttime_ = time;
		ghost_cells_ = ghost_.operator()(tess, cells, time, ts);
		dt_ = dt;
	}

	if (!SmallThermalEnergy(extensive))
		return;
	
	size_t e_index = static_cast<size_t>(entropy_index_);
	assert(entropy_index_<static_cast<int>(extensive.tracers.size())&&entropy_index_>=0);
	
	double new_d=extensive.mass/cd.volumes[index];
	double new_entropy=eos_.dp2s(new_d,NewPressure(extensive, eos_, new_d));
	if(new_entropy*extensive.mass<extensive.tracers[e_index])
		EntropyFix(extensive, cd.volumes[index], eos_, e_index);
	/*else
		return;
	*/
	Vector2D Tgrad = GetTemperatureGrad(index, interp_, cells[index]);
	if ((!BigJump(Tgrad,index,tess,cells,ghost_cells_)&&!NegativeVelocityDivergence(index,interp_,eos_.dp2c(cells[index].density,cells[index].pressure),tess.GetWidth(static_cast<int>(index))))
		||NegativeThermalEnergy(extensive))
	{
		EntropyFix(extensive, cd.volumes[index], eos_, e_index);
		return;
	}
	
/*	Vector2D Tgrad = GetTemperatureGrad(index, interp_, cells[index]);
	if(!PositiveTemperatureDensityGrad(index,interp_,Tgrad))
	{
		EntropyFix(extensive, cd.volumes[index], eos_, e_index);
		return;
	}
	if(!BigJump(Tgrad,index,tess,cells,ghost_cells_))
	{
		EntropyFix(extensive, cd.volumes[index], eos_, e_index);
		return;
	}*/
//	std::cout<<"Have shock"<<std::endl;
}
