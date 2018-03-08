#include "LagrangianExtensiveUpdater3D.hpp"

LagrangianExtensiveUpdater3D::LagrangianExtensiveUpdater3D(LagrangianFlux3D const & lflux, EquationOfState const & eos,
	Ghost3D const & ghost, const vector<pair<const ConditionExtensiveUpdater3D::Condition3D*,
	const ConditionExtensiveUpdater3D::Action3D*>>&sequence) :lflux_(lflux), eos_(eos), ghost_(ghost), sequence_(sequence) {}

void LagrangianExtensiveUpdater3D::operator()(const vector<Conserved3D>& fluxes, const Tessellation3D & tess,
	const double dt, const vector<ComputationalCell3D>& cells, vector<Conserved3D>& extensives, double time,
	TracerStickerNames const & tracerstickernames) const
{
	size_t indexX = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaX")) - tracerstickernames.tracer_names.begin());
	size_t indexY = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaY")) - tracerstickernames.tracer_names.begin());
	size_t indexZ = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaZ")) - tracerstickernames.tracer_names.begin());
	assert(indexX < tracerstickernames.tracer_names.size() && indexY < tracerstickernames.tracer_names.size()
		&& indexZ < tracerstickernames.tracer_names.size());

	size_t N = tess.GetPointNo();
	// Reduce the area tracer
	for (size_t i = 0; i < N; ++i)
	{
		extensives[i].tracers[indexX] *= 0.85;
		extensives[i].tracers[indexY] *= 0.85;
		extensives[i].tracers[indexZ] *= 0.85;
	}

	std::vector<double> oldEk(N, 0), oldEtherm(N, 0), oldE(N, 0);
	for (size_t i = 0; i < N; ++i)
	{
		oldEk[i] = 0.5*ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
		oldEtherm[i] = extensives[i].internal_energy;
		oldE[i] = extensives[i].energy;
	}
	size_t Nfluxes = fluxes.size();
	Conserved3D delta;
	for (size_t i = 0; i < Nfluxes; ++i)
	{
		delta = fluxes[i] * dt*tess.GetArea(i);
		size_t n0 = tess.GetFaceNeighbors(i).first;
		size_t n1 = tess.GetFaceNeighbors(i).second;

		double deltaWs = -lflux_.ws_[i] * tess.GetArea(i) * dt;
		Vector3D normal = normalize(tess.GetMeshPoint(n1) - tess.GetMeshPoint(n0));
		if (lflux_.Lag_calc_[i])
		{
			double p_star = ScalarProd(fluxes[i].momentum, normal);
			double v_star = fluxes[i].energy / p_star;
			double v_new = (v_star - lflux_.ws_[i]);
			if (v_new*v_star > 0)
			{
				if (v_new > 0 && !tess.IsPointOutsideBox(n1) && p_star > 1.2*cells[n1].pressure)
					v_new = std::max(v_new, ScalarProd(cells[n1].velocity, normal));
				if (v_new < 0 && !tess.IsPointOutsideBox(n0) && p_star > 1.2*cells[n0].pressure)
					v_new = std::min(v_new, ScalarProd(cells[n0].velocity, normal));
				delta.energy = p_star*v_new*tess.GetArea(i) * dt;
			}
			else
			{
				if (v_new > 0 && !tess.IsPointOutsideBox(n0))
					delta.energy = tess.GetArea(i) * dt*v_new*cells[n0].pressure;
				else
					if (!tess.IsPointOutsideBox(n1))
						delta.energy = tess.GetArea(i) * dt*v_new*cells[n1].pressure;
					else
						delta.energy = 0;
			}
		}
		if (n0 < N)
		{
			double Ek = 0.5*ScalarProd(extensives[n0].momentum, extensives[n0].momentum) / extensives[n0].mass;
			extensives[n0] -= delta;
			double Eknew = 0.5*ScalarProd(extensives[n0].momentum, extensives[n0].momentum) / extensives[n0].mass;
			extensives[n0].internal_energy -= delta.energy + (Eknew - Ek);
			extensives[n0].tracers[indexX] -= normal.x*deltaWs;
			extensives[n0].tracers[indexY] -= normal.y*deltaWs;
			extensives[n0].tracers[indexZ] -= normal.z*deltaWs;
		}
		if (n1 < N)
		{
			double Ek = 0.5*ScalarProd(extensives[n1].momentum, extensives[n1].momentum) / extensives[n1].mass;
			extensives[n1] += delta;
			double Eknew = 0.5*ScalarProd(extensives[n1].momentum, extensives[n1].momentum) / extensives[n1].mass;
			extensives[n1].internal_energy += delta.energy - (Eknew - Ek);
			extensives[n1].tracers[indexX] -= normal.x*deltaWs;
			extensives[n1].tracers[indexY] -= normal.y*deltaWs;
			extensives[n1].tracers[indexZ] -= normal.z*deltaWs;
		}
	}

	for (size_t i = 0; i < N; ++i)
	{
		double dEtherm = extensives[i].internal_energy - oldEtherm[i];
		double dEk = 0.5*ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass - oldEk[i];
		double dE = extensives[i].energy - oldE[i];
		if (dEtherm*(dE - dEk) > 0)
			if (std::abs(dEtherm) > 0.999 *std::abs(dE - dEk) && std::abs(dEtherm) < 1.001*std::abs(dE - dEk))
				extensives[i].internal_energy = oldEtherm[i] + (dE - dEk);
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
