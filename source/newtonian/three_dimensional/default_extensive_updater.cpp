#include "default_extensive_updater.hpp"
#ifdef RICH_DEBUG
#include <cmath>
#endif

DefaultExtensiveUpdater::DefaultExtensiveUpdater(void){}

void DefaultExtensiveUpdater::operator()(const vector<Conserved3D>& fluxes, const Tessellation3D& tess,
	const double dt, const vector<ComputationalCell3D>& cells, vector<Conserved3D>& extensives, double /*time*/,
	TracerStickerNames const& /*tracerstickernames*/) const
{
	size_t N = tess.GetPointNo();
	size_t Nfluxes = fluxes.size();
	Conserved3D delta;
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
			eo.AddEntry("Face index", i);
			eo.AddEntry("Area", tess.GetArea(i));
			eo.AddEntry("First neigh", n0);
			eo.AddEntry("Second neigh", n1);
			eo.AddEntry("Norg", N);
			eo.AddEntry("Energy flux", fluxes[i].energy);
			eo.AddEntry("Internal Energy flux", fluxes[i].internal_energy);
			eo.AddEntry("Mass flux", fluxes[i].mass);
			eo.AddEntry("Momentum x flux", fluxes[i].momentum.x);
			eo.AddEntry("Momentum y flux", fluxes[i].momentum.y);
			eo.AddEntry("Momentum z flux", fluxes[i].momentum.z);
			for (size_t j = 0; j < delta.tracers.size(); ++j)
				eo.AddEntry("Tracer flux", fluxes[i].tracers[j]);
			eo.AddEntry("Left cell density", cells[n0].density);
			eo.AddEntry("Left cell pressure", cells[n0].pressure);
			eo.AddEntry("Left cell Vx", cells[n0].velocity.x);
			eo.AddEntry("Left cell Vy", cells[n0].velocity.y);
			eo.AddEntry("Left cell Vz", cells[n0].velocity.z);
			eo.AddEntry("Left cell internal energy", cells[n0].internal_energy);
			eo.AddEntry("Left cell ID", cells[n0].ID);
			eo.AddEntry("Right cell density", cells[n1].density);
			eo.AddEntry("Right cell pressure", cells[n1].pressure);
			eo.AddEntry("Right cell Vx", cells[n1].velocity.x);
			eo.AddEntry("Right cell Vy", cells[n1].velocity.y);
			eo.AddEntry("Right cell Vz", cells[n1].velocity.z);
			eo.AddEntry("Right cell internal energy", cells[n1].internal_energy);
			eo.AddEntry("Right cell ID", cells[n1].ID);
			throw eo;
		}
#endif
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
}
