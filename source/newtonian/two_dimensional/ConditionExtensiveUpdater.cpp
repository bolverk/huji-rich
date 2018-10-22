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


