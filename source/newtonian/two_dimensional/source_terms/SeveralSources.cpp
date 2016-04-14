#include "SeveralSources.hpp"

SeveralSources::SeveralSources(vector<SourceTerm*> sources) :
	sources_(sources) {}

SeveralSources::~SeveralSources(void) {}

vector<Extensive> SeveralSources::operator()
(const Tessellation& tess,
	const PhysicalGeometry& pg,
	const CacheData& cd,
	const vector<ComputationalCell>& cells,
	const vector<Extensive>& fluxes,
	const vector<Vector2D>& point_velocities,
	const double t,
	TracerStickerNames const& tracerstickernames) const
{
	vector<Extensive> res(static_cast<size_t>(tess.GetPointNo()));
	size_t ntracer = cells[0].tracers.size();
	for (size_t i = 0; i < res.size(); ++i) 
	{
		res[i].mass = 0;
		res[i].momentum = Vector2D(0, 0);
		res[i].energy = 0;
		res[i].tracers = cells[0].tracers;
		for (size_t j = 0; j < ntracer; ++j)
			res[i].tracers[j] = 0;
	}
	for (size_t i = 0; i < sources_.size(); ++i) 
	{
		const vector<Extensive> diff = (*sources_[i])
			(tess, pg, cd, cells, fluxes, point_velocities, t,tracerstickernames);
		for (size_t j = 0; j < res.size(); ++j) 
		{
			res[j].mass += diff[j].mass;
			res[j].momentum += diff[j].momentum;
			res[j].energy += diff[j].energy;
			for (size_t k = 0; k < ntracer; ++k)
				res[j].tracers[k] += diff[j].tracers[k];
		}
	}
	return res;
}
