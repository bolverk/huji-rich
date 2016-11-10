#include "simple_cell_updater.hpp"
#include "../../misc/lazy_list.hpp"

SimpleCellUpdater::SimpleCellUpdater
(const vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> > sequence) :
	sequence_(sequence), entropy_("Entropy") {}

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
	void regular_update(const EquationOfState& eos, vector<Extensive>& extensives,
		const ComputationalCell& old,
		const CacheData& cd,
		const size_t index,
		ComputationalCell &res,
		size_t entropy_index,
		TracerStickerNames const& tracerstickernames)
	{
		Extensive& extensive = extensives[index];
		const double volume = cd.volumes[index];
		res.density = extensive.mass / volume;
		if (res.density < 0)
			throw UniversalError("Negative density");
		res.velocity = extensive.momentum / extensive.mass;
		const double energy = extensive.energy / extensive.mass -
			0.5*ScalarProd(res.velocity, res.velocity);
		res.stickers = old.stickers;
		for (size_t i = 0; i < extensive.tracers.size(); ++i)
			res.tracers[i] = extensive.tracers[i] / extensive.mass;
		try
		{
			res.pressure = eos.de2p(res.density, energy, res.tracers,tracerstickernames.tracer_names);
			if (entropy_index < res.tracers.size())
			{
				res.tracers[entropy_index] = eos.dp2s(res.density, res.pressure, res.tracers,tracerstickernames.tracer_names);
				extensive.tracers[entropy_index] = res.tracers[entropy_index] * extensive.mass;
			}
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
		TracerStickerNames const & tracerstickernames)
	{
		for (size_t i = 0; i < sequence.size(); ++i)
		{
			if ((*sequence[i].first)(tess, pg, eos, extensives, old, cd, index, tracerstickernames))
			{
				res = (*sequence[i].second)(tess, pg, eos, extensives, old, cd, index, tracerstickernames);
				return;
			}
		}
		regular_update(eos, extensives, old.at(index), cd, index, res, entropyindex,tracerstickernames);
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

	for (size_t i = 0; i < N; ++i)
		update_single(tess, pg, eos, extensives, old, cd, sequence_, i, res[i], tindex,tracerstickernames);
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
