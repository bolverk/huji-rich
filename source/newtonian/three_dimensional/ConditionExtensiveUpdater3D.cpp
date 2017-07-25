#include "ConditionExtensiveUpdater3D.hpp"
#include "../../misc/utils.hpp"
#include <iostream>
#include <cfloat>
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif

ConditionExtensiveUpdater3D::Condition3D::~Condition3D() {}

ConditionExtensiveUpdater3D::Action3D::~Action3D() {}

ConditionExtensiveUpdater3D::~ConditionExtensiveUpdater3D() {}

ConditionExtensiveUpdater3D::ConditionExtensiveUpdater3D(const vector<pair<const Condition3D*, const Action3D*> >& sequence) :
	sequence_(sequence) {}

void ConditionExtensiveUpdater3D::operator()(const vector<Conserved3D>& fluxes,	const Tessellation3D& tess,
	const double dt,const vector<ComputationalCell3D>& cells,vector<Conserved3D>& extensives,double time,
	TracerStickerNames const& tracerstickernames) const
{
	size_t N = tess.GetPointNo();
	size_t Nfluxes = fluxes.size();
	Conserved3D delta;
	for (size_t i = 0; i < Nfluxes; ++i)
	{
		delta = fluxes[i] * dt*tess.GetArea(i);
		size_t n0 = tess.GetFaceNeighbors(i).first;
		size_t n1 = tess.GetFaceNeighbors(i).second;
		if (n0 < N)
		{
			extensives[n0] -= delta;
			extensives[n0].internal_energy -= delta.energy - ScalarProd(cells[n0].velocity, delta.momentum) +
				0.5*ScalarProd(cells[n0].velocity, cells[n0].velocity)*delta.mass;
		}
		if (n1 < N)
		{
			extensives[n1] += delta;
			extensives[n1].internal_energy += delta.energy - ScalarProd(cells[n1].velocity, delta.momentum) +
				0.5*ScalarProd(cells[n1].velocity, cells[n1].velocity)*delta.mass;
		}
	}
	size_t n = tess.GetPointNo();
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < sequence_.size(); ++j)
		{
			if (sequence_[j].first->operator()(i, tess, cells, time, tracerstickernames))
			{
				sequence_[j].second->operator()(fluxes, tess, dt, cells, extensives, i, time, tracerstickernames);
				break;
			}
		}
	}
	extensives.resize(tess.GetPointNo());
}

ColdFlowsUpdate3D::ColdFlowsUpdate3D(EquationOfState const& eos, Ghost3D const& ghost, LinearGauss3D const& interp) :
	eos_(eos), ghost_(ghost), interp_(interp), lasttime_(-DBL_MAX), dt_(0), entropy_index_(1000000), ghost_cells_(boost::container::
		flat_map<size_t, ComputationalCell3D>()) {}

namespace
{
	/*Vector3D GetTemperatureGrad(size_t index, LinearGauss3D const& interp, ComputationalCell3D const& cell)
	{
		Vector3D res;
		Slope3D const& slope = interp.GetSlopesUnlimited()[index];
		res.x = slope.xderivative.pressure / cell.density - slope.xderivative.density*cell.pressure / (cell.density*cell.density);
		res.y = slope.yderivative.pressure / cell.density - slope.yderivative.density*cell.pressure / (cell.density*cell.density);
		res.z = slope.zderivative.pressure / cell.density - slope.zderivative.density*cell.pressure / (cell.density*cell.density);
		return res;
	}

	bool NegativeVelocityDivergence(size_t index, LinearGauss3D const& interp, double cs, double R)
	{
		Slope3D const& slope = interp.GetSlopesUnlimited()[index];
		return (slope.xderivative.velocity.x + slope.yderivative.velocity.y + slope.zderivative.velocity.z)*R < -0.2*cs;
	}
	*/
	/*bool BigJump(Vector3D const& Tgrad, size_t index, Tessellation3D const& tess, vector<ComputationalCell3D> const& cells,
		boost::container::flat_map<size_t, ComputationalCell3D> const& ghost_cells,vector<size_t> &neigh)
	{
		tess.GetNeighbors(index,neigh);
		size_t n = neigh.size();
		size_t N = tess.GetPointNo();
		double maxT = 0, maxP = 0;
		const Vector3D point = tess.GetMeshPoint(index);
		size_t counter = 0;
		for (size_t i = 0; i < n; ++i)
		{
			if (ScalarProd(point - tess.GetMeshPoint(neigh[i]), Tgrad)<0)
			{
				++counter;
				if (neigh[i] < N || !tess.IsPointOutsideBox(neigh[i]))
				{
					maxT = std::max(maxT, cells[neigh[i]].pressure / cells[neigh[i]].density);
					maxP = std::max(maxP, cells[neigh[i]].pressure);
				}
				else
				{
					const boost::container::flat_map<size_t, ComputationalCell3D>::const_iterator it =
						ghost_cells.find(static_cast<size_t>(neigh[i]));
					assert(it != ghost_cells.end());
					maxT = std::max(maxT, it->second.pressure / it->second.density);
					maxP = std::max(maxP, it->second.pressure);
				}
			}
		}
		if (counter == 0)
			return false;
		return (log(maxP) - log(cells[index].pressure))>0.2 ||
			(log(maxT) - log(cells[index].pressure / cells[index].density))>0.1;
	}*/

	/*bool NegativeThermalEnergy(Conserved3D const& cell)
	{
		if (0.50000001*ScalarProd(cell.momentum, cell.momentum) > cell.energy*cell.mass)
			return true;
		else
			return false;
	}
	*/
	bool SmallThermalEnergy(Conserved3D const& cell)
	{
		if (0.51*ScalarProd(cell.momentum, cell.momentum) > cell.energy*cell.mass)
			return true;
		else
			return false;
	}

	bool HighRelativeKineticEnergy(Tessellation3D const& tess, size_t index, vector<Conserved3D> const& cells,
		Conserved3D const& cell)
	{
		vector<size_t> neigh;
		tess.GetNeighbors(index, neigh);
		size_t N = neigh.size();
		double maxDV = 0;
		size_t Norg = tess.GetPointNo();
		Vector3D Vcell = cell.momentum / cell.mass;
		for (size_t i = 0; i < N; ++i)
			if(neigh[i]<Norg||!tess.IsPointOutsideBox(neigh[i]))
				maxDV = std::max(maxDV, abs(Vcell - cells.at(neigh[i]).momentum/cells[neigh[i]].mass));
		return 0.001*maxDV*maxDV*cell.mass > cell.internal_energy;
	}

	void EntropyFix(Conserved3D &extensive, double vol, EquationOfState const& eos, size_t entropy_index)
	{
		const double Ek = 0.5 * ScalarProd(extensive.momentum, extensive.momentum) /  extensive.mass;
		const double density = extensive.mass / vol;
		double Et = eos.dp2e(density, eos.sd2p(extensive.tracers[entropy_index] / extensive.mass, density))
			*extensive.mass;
		if (extensive.internal_energy<0 || Et>extensive.internal_energy*1.03 || Et < extensive.internal_energy*0.97)
		{
			extensive.internal_energy = Et;
			extensive.energy = Ek + Et;
		}
	}

	/*double NewPressure(Conserved3D const& extensive, EquationOfState const& eos, double new_d)
	{
		double new_e = extensive.internal_energy / extensive.mass;
		if (new_e < 0)
			return 0;
		return eos.de2p(new_d, new_e);
	}*/
}

void ColdFlowsUpdate3D::operator()(const vector<Conserved3D>& /*fluxes*/, const Tessellation3D& tess, const double dt,
	const vector<ComputationalCell3D>& cells, vector<Conserved3D> &extensives, size_t index, double time, TracerStickerNames const& ts)
	const
{
	vector<size_t> neigh;
	if (entropy_index_ > 100000)
	{
		entropy_index_ = static_cast<int>(lower_bound(ts.tracer_names.begin(), ts.tracer_names.end(), string("Entropy")) - ts.tracer_names.begin());
		lasttime_ = time;
		ghost_.operator()(tess, cells, time, ts, ghost_cells_);
#ifdef RICH_MPI
		MPI_exchange_data(tess, extensives, true);
#endif
		dt_ = dt;
	}
	if (lasttime_ < time || dt < dt_ || dt > dt_)
	{
		lasttime_ = time;
		ghost_.operator()(tess, cells, time, ts, ghost_cells_);
		dt_ = dt;
#ifdef RICH_MPI
		MPI_exchange_data(tess, extensives, true);
#endif
	}

	if (!SmallThermalEnergy(extensives[index]))
		return;

	assert(entropy_index_ < extensives[index].tracers.size());

	double new_d = extensives[index].mass / tess.GetVolume(index);
	assert(new_d > 0);
	if (!HighRelativeKineticEnergy(tess, index, extensives,extensives[index]))
		return;
	else
	{
		EntropyFix(extensives[index], tess.GetVolume(index), eos_, entropy_index_);
		return;
	}
	/*
	double new_entropy = eos_.dp2s(new_d, NewPressure(extensive, eos_, new_d));
	if (new_entropy*extensive.mass < extensive.tracers[entropy_index_])
	{
		EntropyFix(extensive, tess.GetVolume(index), eos_, entropy_index_);
		return;
	}

	Vector3D Tgrad = GetTemperatureGrad(index, interp_, cells[index]);
	if ((!BigJump(Tgrad, index, tess, cells, ghost_cells_, neigh) &&
		!NegativeVelocityDivergence(index, interp_, eos_.dp2c(cells[index].density, cells[index].pressure), tess.GetWidth(index)))
		|| NegativeThermalEnergy(extensive))
	{
		EntropyFix(extensive, tess.GetVolume(index), eos_, entropy_index_);
		return;
	}
	return;*/
}

void RegularExtensiveUpdate3D::operator()(const vector<Conserved3D>& /*fluxes*/, const Tessellation3D& /*tess*/, const double /*dt*/,
	const vector<ComputationalCell3D>& /*cells*/, vector<Conserved3D> &/*extensives*/, size_t /*index*/, double /*time*/,
	TracerStickerNames const& /*tracerstickernames*/)const
{
	return;
}