#include "default_extensive_updater.hpp"
#ifdef RICH_DEBUG
#include <cmath>
#endif

DefaultExtensiveUpdater::DefaultExtensiveUpdater(void){}

void DefaultExtensiveUpdater::operator()(const vector<Conserved3D>& fluxes, const Tessellation3D& tess,
	const double dt, const vector<ComputationalCell3D>& cells, vector<Conserved3D>& extensives, double /*time*/,
					 TracerStickerNames const& /*tracerstickernames*/, std::vector<Vector3D> const& /*face_vel*/,
					 std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > const& /*interp_values*/) const
{
	size_t N = tess.GetPointNo();
	size_t Nfluxes = fluxes.size();
	Conserved3D delta;
	std::vector<double> oldEk(N, 0), oldEtherm(N, 0), oldE(N, 0);
	for (size_t i = 0; i < N; ++i)
	{
		oldEk[i] = 0.5*ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
		oldEtherm[i] = extensives[i].internal_energy;
		oldE[i] = extensives[i].energy;
	}
	for (size_t i = 0; i < Nfluxes; ++i)
	{
		delta = fluxes[i] * dt*tess.GetArea(i);
		size_t n0 = tess.GetFaceNeighbors(i).first;
		size_t n1 = tess.GetFaceNeighbors(i).second;
#ifdef RICH_DEBUG
		bool good = true;
		if (!std::isfinite(fluxes[i].energy))
			good = false;
		if (!std::isfinite(fluxes[i].internal_energy))
			good = false;
		if (!std::isfinite(fluxes[i].momentum.x))
			good = false;
		if (!std::isfinite(fluxes[i].momentum.y))
			good = false;
		if (!std::isfinite(fluxes[i].momentum.z))
			good = false;
		if (!std::isfinite(fluxes[i].mass))
			good = false;
		for (size_t j = 0; j < delta.tracers.size(); ++j)
		{
			if (!std::isfinite(fluxes[i].tracers[j]))
				good = false;
		}
		if (!good)
		{
			UniversalError eo("Bad flux");
			eo.addEntry("Face index", i);
			eo.addEntry("Area", tess.GetArea(i));
			eo.addEntry("First neigh", n0);
			eo.addEntry("Second neigh", n1);
			eo.addEntry("Norg", N);
			eo.addEntry("Energy flux", fluxes[i].energy);
			eo.addEntry("Internal Energy flux", fluxes[i].internal_energy);
			eo.addEntry("Mass flux", fluxes[i].mass);
			eo.addEntry("Momentum x flux", fluxes[i].momentum.x);
			eo.addEntry("Momentum y flux", fluxes[i].momentum.y);
			eo.addEntry("Momentum z flux", fluxes[i].momentum.z);
			for (size_t j = 0; j < delta.tracers.size(); ++j)
				eo.addEntry("Tracer flux", fluxes[i].tracers[j]);
			eo.addEntry("Left cell density", cells[n0].density);
			eo.addEntry("Left cell pressure", cells[n0].pressure);
			eo.addEntry("Left cell Vx", cells[n0].velocity.x);
			eo.addEntry("Left cell Vy", cells[n0].velocity.y);
			eo.addEntry("Left cell Vz", cells[n0].velocity.z);
			eo.addEntry("Left cell internal energy", cells[n0].internal_energy);
			eo.addEntry("Left cell ID", cells[n0].ID);
			eo.addEntry("Right cell density", cells[n1].density);
			eo.addEntry("Right cell pressure", cells[n1].pressure);
			eo.addEntry("Right cell Vx", cells[n1].velocity.x);
			eo.addEntry("Right cell Vy", cells[n1].velocity.y);
			eo.addEntry("Right cell Vz", cells[n1].velocity.z);
			eo.addEntry("Right cell internal energy", cells[n1].internal_energy);
			eo.addEntry("Right cell ID", cells[n1].ID);
			throw eo;
		}
#endif
		delta.internal_energy = 0;
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
	for (size_t i = 0; i < N; ++i)
	{
		double dEtherm = extensives[i].internal_energy - oldEtherm[i];
		double Eknew = 0.5*ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
		double dEk = Eknew - oldEk[i];
		double dE = extensives[i].energy - oldE[i];
		if (dEtherm*(dE - dEk) > 0)
		{
			if (std::abs(dEtherm) > 0.95 *std::abs(dE - dEk) && std::abs(dEtherm) < 1.05*std::abs(dE - dEk))
				extensives[i].internal_energy = extensives[i].energy - Eknew;
			else
				extensives[i].energy = extensives[i].internal_energy + Eknew;
		}
		else
			extensives[i].energy = extensives[i].internal_energy + Eknew;
	}
	extensives.resize(tess.GetPointNo());
}
