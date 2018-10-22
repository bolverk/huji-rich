#include "default_cell_updater.hpp"
#include "../../misc/utils.hpp"
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif

namespace
{
	void EntropyFix(EquationOfState const& eos, ComputationalCell3D &res, size_t entropy_index, TracerStickerNames const& tracerstickernames, double &energy,
		Conserved3D &extensive)
	{
		double new_pressure = eos.sd2p(res.tracers[entropy_index], res.density, res.tracers, tracerstickernames.tracer_names);
		res.pressure = new_pressure;
		double de = eos.dp2e(res.density, res.pressure) - energy;
		energy += de;
		extensive.energy += de * extensive.mass;
		extensive.internal_energy += de * extensive.mass;
	}

	bool HighRelativeKineticEnergy(Tessellation3D const& tess, size_t index, vector<Conserved3D> const& cells,
		Conserved3D const& cell)
	{
		std::vector<size_t> neigh;
		tess.GetNeighbors(index, neigh);
		size_t N = neigh.size();
		double maxDV = 0;
		size_t Norg = tess.GetPointNo();
		Vector3D Vcell = cell.momentum / cell.mass;
		double Et = cell.energy / cell.mass - 0.5*ScalarProd(Vcell, Vcell);
		for (size_t i = 0; i < N; ++i)
			if (neigh[i] < Norg || !tess.IsPointOutsideBox(neigh[i]))
				maxDV = std::max(maxDV, abs(Vcell - cells.at(neigh[i]).momentum / cells[neigh[i]].mass));
		return 0.005*maxDV*maxDV > Et;
	}

}

DefaultCellUpdater::DefaultCellUpdater(void):entropy_index_(9999999){}

void DefaultCellUpdater::operator()(vector<ComputationalCell3D> &res, EquationOfState const& eos,
	const Tessellation3D& tess, vector<Conserved3D>& extensives, 
	TracerStickerNames const& tracerstickernames) const
{
	entropy_index_ = tracerstickernames.tracer_names.size();
	vector<string>::const_iterator it = binary_find(tracerstickernames.tracer_names.begin(),
		tracerstickernames.tracer_names.end(), string("Entropy"));
	if (it != tracerstickernames.tracer_names.end())
		entropy_index_ = static_cast<size_t>(it - tracerstickernames.tracer_names.begin());
#ifdef RICH_MPI
	if (entropy_index_ < tracerstickernames.tracer_names.size())
		MPI_exchange_data(tess, extensives, true);
#endif
	size_t Nloop = tess.GetPointNo();
	size_t Ntracers = res.at(0).tracers.size();
	res.resize(Nloop);
	for (size_t i = 0; i < Nloop; ++i)
	{
		Conserved3D& extensive = extensives[i];
		const double vol = tess.GetVolume(i);
		res[i].density = extensive.mass / vol;
		res[i].velocity = extensive.momentum / extensive.mass;
		double energy = extensive.internal_energy / extensive.mass;
		if (energy < 0)
		{
			energy = extensive.energy / extensive.mass - 0.5*ScalarProd(res[i].velocity, res[i].velocity);
			extensive.internal_energy = energy*extensive.mass;
		}
		extensive.energy = extensive.mass*(energy + 0.5*ScalarProd(res[i].velocity, res[i].velocity));
		for (size_t j = 0; j < Ntracers; ++j)
			res[i].tracers[j] = extensive.tracers[j] / extensive.mass;
		try
		{
			// Entropy fix if needed
			if (entropy_index_ < tracerstickernames.tracer_names.size())
			{
				// Do we have a negative thermal energy?
				if (energy < 0)
				{
					EntropyFix(eos, res[i], entropy_index_, tracerstickernames, energy, extensive);
				}
				else
				{
					double new_pressure = eos.de2p(res[i].density, energy, res[i].tracers, tracerstickernames.tracer_names);
					double new_entropy = eos.dp2s(res[i].density, new_pressure, res[i].tracers, tracerstickernames.tracer_names);
					// Did we lose entropy? If yes then correct for it since it is unphysical
					if (new_entropy < res[i].tracers[entropy_index_])
					{
						EntropyFix(eos, res[i], entropy_index_, tracerstickernames, energy, extensive);
					}
					else
					{
						// Is the kinetic energy small?
						if (energy*extensive.mass < 0.005*extensive.energy)
						{
							if (!HighRelativeKineticEnergy(tess, i, extensives, extensive))
							{
								EntropyFix(eos, res[i], entropy_index_, tracerstickernames, energy, extensive);
							}
						}
						else
						{
							// We don't need the entropy fix, update entropy
							res[i].pressure = new_pressure;
							res[i].tracers[entropy_index_] = new_entropy;
							extensive.tracers[entropy_index_] = new_entropy * extensive.mass;
						}
					}
				}
			}
			else
			{
				res[i].pressure = eos.de2p(res[i].density, energy, res[i].tracers, tracerstickernames.tracer_names);
				res[i].internal_energy = energy;
			}
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Cell index", static_cast<double>(i));
			eo.AddEntry("Cell mass", extensive.mass);
			eo.AddEntry("Cell x momentum", extensive.momentum.x);
			eo.AddEntry("Cell y momentum", extensive.momentum.y);
			eo.AddEntry("Cell zy momentum", extensive.momentum.z);
			eo.AddEntry("Cell energy", extensive.energy);
			throw;
		}
		if (res[i].density < 0 || energy < 0)
		{
			UniversalError eo("Negative quantity in cell update");
			eo.AddEntry("Cell index", static_cast<double>(i));
			eo.AddEntry("Cell mass", extensive.mass);
			eo.AddEntry("Cell x momentum", extensive.momentum.x);
			eo.AddEntry("Cell y momentum", extensive.momentum.y);
			eo.AddEntry("Cell z momentum", extensive.momentum.z);
			eo.AddEntry("Cell x location", tess.GetMeshPoint(i).x);
			eo.AddEntry("Cell y location", tess.GetMeshPoint(i).y);
			eo.AddEntry("Cell z location", tess.GetMeshPoint(i).z);
			eo.AddEntry("Cell volume", vol);
			eo.AddEntry("Cell energy", extensive.energy);
			eo.AddEntry("Cell thermal energy per unit mass", energy);
			throw eo;
		}
	}
}
