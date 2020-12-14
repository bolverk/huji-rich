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
		res.internal_energy = energy;
		extensive.energy += de * extensive.mass;
		extensive.internal_energy += de * extensive.mass;
	}

	void EntropyFixSR(EquationOfState const& eos, ComputationalCell3D &res, size_t entropy_index, TracerStickerNames const& tracerstickernames, Conserved3D &extensive, double vol)
	{
		double new_pressure = eos.sd2p(res.tracers[entropy_index], res.density, res.tracers, tracerstickernames.tracer_names);
		res.pressure = new_pressure;
		double energy = eos.dp2e(res.density, res.pressure);
		double gamma = std::sqrt(1.0 / (1 - ScalarProd(res.velocity, res.velocity)));
		extensive.energy = extensive.mass*(energy*gamma + gamma - 1) - res.pressure*vol;
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

	void regular_update(std::vector<ComputationalCell3D> &res, std::vector<Conserved3D> & extensives,
		Tessellation3D const& tess, size_t entropy_index, TracerStickerNames const& tsn,
		EquationOfState const& eos)
	{
		size_t Nloop = tess.GetPointNo();
		size_t Ntracers = tsn.tracer_names.size();
		for (size_t i = 0; i < Nloop; ++i)
		{
			try
			{
				Conserved3D& extensive = extensives[i];
				const double vol = tess.GetVolume(i);
				res[i].density = extensive.mass / vol;
				res[i].velocity = extensive.momentum / extensive.mass;
				double energy = extensive.internal_energy / extensive.mass;
				extensive.energy = extensive.mass*(energy + 0.5*ScalarProd(res[i].velocity, res[i].velocity));
				for (size_t j = 0; j < Ntracers; ++j)
					res[i].tracers[j] = extensive.tracers[j] / extensive.mass;
				// Entropy fix if needed
				if (entropy_index < tsn.tracer_names.size())
				{
					try
					{
						// Do we have a negative thermal energy?
						if (energy < 0)
						{
							EntropyFix(eos, res[i], entropy_index, tsn, energy, extensive);
						}
						else
						{
							// Is the kinetic energy small?
							if ((energy*extensive.mass < 0.005*extensive.energy) &&
								HighRelativeKineticEnergy(tess, i, extensives, extensive))
							{
								EntropyFix(eos, res[i], entropy_index, tsn, energy, extensive);
							}
							else
							{
								double new_pressure = eos.de2p(res[i].density, energy, res[i].tracers, tsn.tracer_names);
								double new_entropy = eos.dp2s(res[i].density, new_pressure, res[i].tracers, tsn.tracer_names);
								// We don't need the entropy fix, update entropy
								res[i].internal_energy = energy;
								res[i].pressure = new_pressure;
								res[i].tracers[entropy_index] = new_entropy;
								extensive.tracers[entropy_index] = new_entropy * extensive.mass;
							}
						}
					}
					catch (UniversalError &eo)
					{
						eo.AddEntry("Cell index", static_cast<double>(i));
						eo.AddEntry("Cell mass", extensives[i].mass);
						eo.AddEntry("Cell x momentum", extensives[i].momentum.x);
						eo.AddEntry("Cell y momentum", extensives[i].momentum.y);
						eo.AddEntry("Cell z momentum", extensives[i].momentum.z);
						eo.AddEntry("Cell x location", tess.GetMeshPoint(i).x);
						eo.AddEntry("Cell y location", tess.GetMeshPoint(i).y);
						eo.AddEntry("Cell z location", tess.GetMeshPoint(i).z);
						eo.AddEntry("Cell volume", vol);
						eo.AddEntry("Cell energy", extensives[i].energy);
						eo.AddEntry("Cell thermal energy per unit mass", energy);
						eo.AddEntry("Cell id", static_cast<double>(res[i].ID));
						throw eo;
					}
				}
				else
				{
					res[i].pressure = eos.de2p(res[i].density, energy, res[i].tracers, tsn.tracer_names);
					res[i].internal_energy = energy;
				}
				if (!(res[i].density > 0) || !(res[i].pressure > 0) || (!std::isfinite(fastabs(extensives[i].momentum))))
				{
					UniversalError eo("Negative quantity in cell update");
					eo.AddEntry("Cell index", static_cast<double>(i));
					eo.AddEntry("Cell mass", extensives[i].mass);
					eo.AddEntry("Cell x momentum", extensives[i].momentum.x);
					eo.AddEntry("Cell y momentum", extensives[i].momentum.y);
					eo.AddEntry("Cell z momentum", extensives[i].momentum.z);
					eo.AddEntry("Cell x location", tess.GetMeshPoint(i).x);
					eo.AddEntry("Cell y location", tess.GetMeshPoint(i).y);
					eo.AddEntry("Cell z location", tess.GetMeshPoint(i).z);
					eo.AddEntry("Cell volume", vol);
					eo.AddEntry("Cell energy", extensives[i].energy);
					eo.AddEntry("Cell thermal energy per unit mass", energy);
					eo.AddEntry("Cell id", static_cast<double>(res[i].ID));
					throw eo;
				}
			}
			catch (UniversalError &eo)
			{
				eo.AddEntry("Cell index", static_cast<double>(i));
				eo.AddEntry("Cell mass", extensives[i].mass);
				eo.AddEntry("Cell x momentum", extensives[i].momentum.x);
				eo.AddEntry("Cell y momentum", extensives[i].momentum.y);
				eo.AddEntry("Cell z momentum", extensives[i].momentum.z);
				eo.AddEntry("Cell energy", extensives[i].energy);
				throw eo;
			}
		}
	}

	void regular_updateSR(std::vector<ComputationalCell3D> &res, std::vector<Conserved3D> & extensives,
		Tessellation3D const& tess, size_t entropy_index, TracerStickerNames const& tsn,
		EquationOfState const& eos, double G)
	{
		size_t Nloop = tess.GetPointNo();
		size_t Ntracers = tsn.tracer_names.size();
		for (size_t i = 0; i < Nloop; ++i)
		{
			try
			{
				double v = GetVelocity(extensives[i], G);
				const double volume = 1.0 / tess.GetVolume(i);
				if (fastabs(extensives[i].momentum)*1e8 < extensives[i].mass)
				{
					res[i].velocity = extensives[i].momentum / extensives[i].mass;
					v = abs(res[i].velocity);
				}
				else
					res[i].velocity = v * extensives[i].momentum / abs(extensives[i].momentum);
				double gamma_1 = std::sqrt(1 - ScalarProd(res[i].velocity, res[i].velocity));
				res[i].density = extensives[i].mass *gamma_1*volume;
				//				double r = fastabs(tess.GetMeshPoint(i));
				if (res[i].density < 0)
					throw UniversalError("Negative density");
				for (size_t j = 0; j < extensives[i].tracers.size(); ++j)
					res[i].tracers[j] = extensives[i].tracers[j] / extensives[i].mass;
				if (fastabs(res[i].velocity) < 1e-5)
					res[i].pressure = (G - 1)*((extensives[i].energy - ScalarProd(extensives[i].momentum, res[i].velocity))*volume
						+ (0.5*ScalarProd(res[i].velocity, res[i].velocity))*res[i].density);
				else
					res[i].pressure = (G - 1)*(extensives[i].energy*volume - ScalarProd(extensives[i].momentum, res[i].velocity)*volume
						+ (1.0 / gamma_1 - 1)*res[i].density);
				// Entropy fix if needed
				if (entropy_index < Ntracers)
				{
					if (!(extensives[i].tracers[entropy_index] > 0))
					{
						UniversalError eo("Negative entropy");
						eo.AddEntry("Extensive Entropy", extensives[i].tracers[entropy_index]);
						throw eo;
					}
					// Do we have a negative thermal energy?
					if (res[i].pressure < 0)
					{
						EntropyFixSR(eos, res[i], entropy_index, tsn, extensives[i], 1.0 / volume);
					}
					else
					{
						double new_entropy = eos.dp2s(res[i].density, res[i].pressure, res[i].tracers, tsn.tracer_names);
						// We don't need the entropy fix, update entropy
						res[i].tracers[entropy_index] = new_entropy;
						extensives[i].tracers[entropy_index] = new_entropy * extensives[i].mass;
					}
				}
				if (!(res[i].density > 0) || !(res[i].pressure > 0)||(!std::isfinite(fastabs(extensives[i].momentum))))
				{
					UniversalError eo("Negative quantity in cell update");
					throw eo;
				}
				res[i].internal_energy = eos.dp2e(res[i].density, res[i].pressure,res[i].tracers,tsn.tracer_names);
			}
			catch (UniversalError &eo)
			{
			  eo.AddEntry("Cell ID", static_cast<double>(res[i].ID));
				eo.AddEntry("Cell volume", tess.GetVolume(i));
				eo.AddEntry("Cell index", static_cast<double>(i));
				eo.AddEntry("Cell mass", extensives[i].mass);
				eo.AddEntry("Cell x momentum", extensives[i].momentum.x);
				eo.AddEntry("Cell y momentum", extensives[i].momentum.y);
				eo.AddEntry("Cell z momentum", extensives[i].momentum.z);
				eo.AddEntry("Cell x location", tess.GetMeshPoint(i).x);
				eo.AddEntry("Cell y location", tess.GetMeshPoint(i).y);
				eo.AddEntry("Cell z location", tess.GetMeshPoint(i).z);
				eo.AddEntry("Cell energy", extensives[i].energy);
				eo.AddEntry("Number of tracers", static_cast<double>(Ntracers));
				for (size_t j = 0; j < Ntracers; ++j)
					eo.AddEntry("Tracer extensive value", extensives[i].tracers[j]);
				eo.AddEntry("Entropy index", static_cast<double>(entropy_index));
				throw;
			}
		}
	}

}

DefaultCellUpdater::DefaultCellUpdater(bool SR, double G) :SR_(SR), G_(G), entropy_index_(9999999) {}

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
	{
		Conserved3D edummy;
		MPI_exchange_data(tess, extensives, true,&edummy);
	}
#endif
	if (!SR_)
		regular_update(res, extensives, tess, entropy_index_, tracerstickernames, eos);
	else
		regular_updateSR(res, extensives, tess, entropy_index_, tracerstickernames, eos, G_);
#ifdef RICH_MPI
	if (entropy_index_ < tracerstickernames.tracer_names.size())
		extensives.resize(tess.GetPointNo());
#endif

}
