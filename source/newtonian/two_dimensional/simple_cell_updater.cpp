#include "simple_cell_updater.hpp"
#include "../../misc/lazy_list.hpp"
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif

SimpleCellUpdater::SimpleCellUpdater
(const vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> > sequence, bool SR, double G) :
	sequence_(sequence), SR_(SR), G_(G), entropy_("Entropy") {}

SimpleCellUpdater::~SimpleCellUpdater(void)
{
	for (size_t i = 0; i < sequence_.size(); ++i) {
		delete sequence_[i].first;
		delete sequence_[i].second;
	}
}

SimpleCellUpdater::Condition::~Condition(void) {}

SimpleCellUpdater::Action::~Action(void) {}

namespace
{
	void EntropyFix(EquationOfState const& eos, ComputationalCell &res, size_t entropy_index, TracerStickerNames const& tracerstickernames, double &energy,
		Extensive &extensive)
	{
		double new_pressure = eos.sd2p(res.tracers[entropy_index], res.density, res.tracers, tracerstickernames.tracer_names);
		res.pressure = new_pressure;
		double de = eos.dp2e(res.density, res.pressure) - energy;
		energy += de;
		extensive.energy += de * extensive.mass;
	}

	void EntropyFixSR(EquationOfState const& eos, ComputationalCell &res, size_t entropy_index, TracerStickerNames const& tracerstickernames, double &energy,
		Extensive &extensive, double vol)
	{
		double new_pressure = eos.sd2p(res.tracers[entropy_index], res.density, res.tracers, tracerstickernames.tracer_names);
		res.pressure = new_pressure;
		energy = eos.dp2e(res.density, res.pressure);
		double gamma = std::sqrt(1.0 / (1 - ScalarProd(res.velocity, res.velocity)));
		extensive.energy = extensive.mass*(energy*gamma + gamma - 1) - res.pressure*vol;
	}

	bool HighRelativeKineticEnergy(Tessellation const& tess, size_t index, vector<Extensive> const& cells,
		Extensive const& cell)
	{
		std::vector<int> neigh;
		tess.GetNeighbors(static_cast<int>(index), neigh);
		size_t N = neigh.size();
		double maxDV = 0;
		size_t Norg = static_cast<size_t>(tess.GetPointNo());
		Vector2D Vcell = cell.momentum / cell.mass;
		double Et = cell.energy / cell.mass - 0.5*ScalarProd(Vcell, Vcell);
		for (size_t i = 0; i < N; ++i)
			if (neigh[i] < static_cast<int>(Norg) || tess.GetOriginalIndex(static_cast<int>(neigh[i]))>static_cast<int>(Norg))
				maxDV = std::max(maxDV, abs(Vcell - cells.at(static_cast<size_t>(neigh[i])).momentum / cells[static_cast<size_t>(neigh[i])].mass));
		return 0.005*maxDV*maxDV > Et;
	}

	void regular_update(const EquationOfState& eos, vector<Extensive>& extensives,
		const ComputationalCell& old,
		const CacheData& cd,
		const size_t index,
		ComputationalCell &res,
		size_t entropy_index,
		TracerStickerNames const& tracerstickernames,
		Tessellation const& tess)
	{
		Extensive& extensive = extensives[index];
		const double volume = cd.volumes[index];
		res.density = extensive.mass / volume;
		if (res.density < 0)
			throw UniversalError("Negative density");
		res.velocity = extensive.momentum / extensive.mass;
		double energy = extensive.energy / extensive.mass -
			0.5*ScalarProd(res.velocity, res.velocity);
		res.stickers = old.stickers;
		for (size_t i = 0; i < extensive.tracers.size(); ++i)
			res.tracers[i] = extensive.tracers[i] / extensive.mass;
		try
		{
			// Entropy fix if needed
			if (entropy_index < res.tracers.size())
			{
				// Do we have a negative thermal energy?
				if (energy < 0)
				{
					EntropyFix(eos, res, entropy_index, tracerstickernames, energy, extensive);
				}
				else
				{
					// Is the kinetic energy small?
					if ((energy*extensive.mass < 0.005*extensive.energy) &&
						!HighRelativeKineticEnergy(tess, index, extensives, extensive))
					{
						EntropyFix(eos, res, entropy_index, tracerstickernames, energy, extensive);
					}
					else
					{
						double new_pressure = eos.de2p(res.density, energy, res.tracers, tracerstickernames.tracer_names);
						double new_entropy = eos.dp2s(res.density, new_pressure, res.tracers, tracerstickernames.tracer_names);
						// We don't need the entropy fix, update entropy
						res.pressure = new_pressure;
						res.tracers[entropy_index] = new_entropy;
						extensive.tracers[entropy_index] = new_entropy * extensive.mass;
					}
				}
			}
			else
				res.pressure = eos.de2p(res.density, energy, res.tracers, tracerstickernames.tracer_names);
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Cell index", static_cast<double>(index));
			eo.AddEntry("Cell mass", extensive.mass);
			eo.AddEntry("Cell x momentum", extensive.momentum.x);
			eo.AddEntry("Cell y momentum", extensive.momentum.y);
			eo.AddEntry("Cell energy", extensive.energy);
			throw;
		}
	}

	void regular_updateSR(const EquationOfState& eos, vector<Extensive>& extensives,
		const ComputationalCell& old,
		const CacheData& cd,
		const size_t index,
		ComputationalCell &res,
		size_t entropy_index,
		TracerStickerNames const& tracerstickernames,
		Tessellation const& /*tess*/, double G)
	{
		Extensive& extensive = extensives[index];
		double v = GetVelocity(extensive, G);
		const double volume = 1.0 / cd.volumes[index];
		res.velocity = (fastabs(extensive.momentum)*1e8 < extensive.mass) ? extensive.momentum / extensive.mass : v * extensive.momentum / abs(extensive.momentum);
		double gamma_1 = std::sqrt(1 - ScalarProd(res.velocity, res.velocity));
		res.density = extensive.mass *gamma_1*volume;
		if (res.density < 0)
			throw UniversalError("Negative density");
		res.stickers = old.stickers;
		for (size_t i = 0; i < extensive.tracers.size(); ++i)
			res.tracers[i] = extensive.tracers[i] / extensive.mass;
		if (fastabs(res.velocity) < 1e-5)
			res.pressure = (G - 1)*((extensive.energy - ScalarProd(extensive.momentum, res.velocity))*volume
				+ (0.5*ScalarProd(res.velocity, res.velocity))*res.density);
		else
			res.pressure = (G - 1)*(extensive.energy*volume - ScalarProd(extensive.momentum, res.velocity)*volume
				+ (1.0 / gamma_1 - 1)*res.density);
		// Entropy fix if needed
		if (entropy_index < res.tracers.size())
		{
			// Do we have a negative thermal energy?
			if (res.pressure < 0)
			{
				double enthalpy = 0;
				EntropyFixSR(eos, res, entropy_index, tracerstickernames, enthalpy, extensive, cd.volumes[index]);
			}
			else
			{
				double new_entropy = eos.dp2s(res.density, res.pressure, res.tracers, tracerstickernames.tracer_names);
				// We don't need the entropy fix, update entropy
				res.tracers[entropy_index] = new_entropy;
				extensive.tracers[entropy_index] = new_entropy * extensive.mass;
			}
		}
	}

	void update_single(const Tessellation& tess,
		const PhysicalGeometry& pg,
		const EquationOfState& eos,
		vector<Extensive>& extensives,
		const vector<ComputationalCell>& old,
		const CacheData& cd,
		const vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> >& sequence,
		const size_t index,
		ComputationalCell &res,
		size_t entropyindex,
		TracerStickerNames const & tracerstickernames, bool SR, double G)
	{
		for (size_t i = 0; i < sequence.size(); ++i)
		{
			if ((*sequence[i].first)(tess, pg, eos, extensives, old, cd, index, tracerstickernames))
			{
				res = (*sequence[i].second)(tess, pg, eos, extensives, old, cd, index, tracerstickernames);
				return;
			}
		}
		if (!SR)
			regular_update(eos, extensives, old.at(index), cd, index, res, entropyindex, tracerstickernames, tess);
		else
			regular_updateSR(eos, extensives, old.at(index), cd, index, res, entropyindex, tracerstickernames, tess, G);
	}
}

vector<ComputationalCell> SimpleCellUpdater::operator()
(const Tessellation& tess,
	const PhysicalGeometry& pg,
	const EquationOfState& eos,
	vector<Extensive>& extensives,
	const vector<ComputationalCell>& old,
	const CacheData& cd,
	TracerStickerNames const& tracerstickernames) const
{
	size_t N = static_cast<size_t>(tess.GetPointNo());
	vector<ComputationalCell> res(N, old[0]);

	size_t tindex = old[0].tracers.size();
	vector<string>::const_iterator it = binary_find(tracerstickernames.tracer_names.begin(),
		tracerstickernames.tracer_names.end(), entropy_);
	if (it != tracerstickernames.tracer_names.end())
		tindex = static_cast<size_t>(it - tracerstickernames.tracer_names.begin());
#ifdef RICH_MPI
	if (tindex < old[0].tracers.size())
		MPI_exchange_data(tess, extensives, true);
#endif
	for (size_t i = 0; i < N; ++i)
		update_single(tess, pg, eos, extensives, old, cd, sequence_, i, res[i], tindex, tracerstickernames, SR_, G_);
#ifdef RICH_MPI
	if (tindex < old[0].tracers.size())
		extensives.resize(static_cast<size_t>(tess.GetPointNo()));
#endif
	return res;
}

HasSticker::HasSticker
(const string& sticker_name) :
	sticker_name_(sticker_name) {}

bool HasSticker::operator()
(const Tessellation& /*tess*/,
	const PhysicalGeometry& /*pg*/,
	const EquationOfState& /*eos*/,
	const vector<Extensive>& /*extensives*/,
	const vector<ComputationalCell>& cells,
	const CacheData& /*cd*/,
	const size_t index,
	TracerStickerNames const& tracerstickernames) const
{
	vector<string>::const_iterator it = binary_find(tracerstickernames.sticker_names.begin(), tracerstickernames.sticker_names.end(),
		sticker_name_);
	assert(it != tracerstickernames.sticker_names.end());
	return cells[index].stickers[static_cast<size_t>(it - tracerstickernames.sticker_names.begin())];
}

SkipUpdate::SkipUpdate(void) {}

ComputationalCell SkipUpdate::operator()
(const Tessellation& /*tess*/,
	const PhysicalGeometry& /*pg*/,
	const EquationOfState& /*eos*/,
	const vector<Extensive>& /*extensives*/,
	const vector<ComputationalCell>& cells,
	const CacheData& /*cd*/,
	const size_t index,
	TracerStickerNames const& /*tracerstickernames*/) const
{
	return cells[index];
}
