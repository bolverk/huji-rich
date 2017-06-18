#include "PCM.hpp"

PCM::PCM(GhostPointGenerator const & ghost) :ghost_(ghost){}

void PCM::operator()(const Tessellation & tess, const vector<ComputationalCell>& cells, double time,
	vector<pair<ComputationalCell, ComputationalCell> >& res, TracerStickerNames const & tracerstickersnames,
	CacheData const& /*cd*/) const
{
	res.resize(static_cast<size_t>(tess.GetTotalSidesNumber()));
	boost::container::flat_map<size_t, ComputationalCell> ghost_cells = ghost_.operator()(tess,
		cells, time, tracerstickersnames);
	int Npoints = tess.GetPointNo();
	size_t Nedges = res.size();
	for (size_t i = 0; i < Nedges; ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(i));
		if (edge.neighbors.first < Npoints)
			res[i].first = cells[static_cast<size_t>(edge.neighbors.first)];
		else
#ifdef RICH_MPI
			if (tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first))
				res[i].first = cells.at(static_cast<size_t>(edge.neighbors.first));
			else
#endif
				res[i].first = ghost_cells[static_cast<size_t>(edge.neighbors.first)];
		if (edge.neighbors.second < Npoints)
			res[i].second = cells[static_cast<size_t>(edge.neighbors.second)];
		else
#ifdef RICH_MPI
			if (tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first))
				res[i].second = cells.at(static_cast<size_t>(edge.neighbors.second));
			else
#endif
				res[i].second = ghost_cells[static_cast<size_t>(edge.neighbors.second)];

	}
}
