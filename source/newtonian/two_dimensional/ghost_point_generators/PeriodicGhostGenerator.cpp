#include "PeriodicGhostGenerator.hpp"

boost::container::flat_map<size_t, ComputationalCell> PeriodicGhostGenerator::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells, double /*time*/, TracerStickerNames const& /*tracerstickernames*/) const
{
	boost::container::flat_map<size_t, ComputationalCell> res;
	vector<std::pair<size_t, size_t> > ghosts = GetOuterEdgesIndeces(tess);
	for (size_t i = 0; i < ghosts.size(); ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(ghosts[i].first));
		ComputationalCell ctemp = cells[static_cast<size_t>(ghosts[i].second == 2 ? 
			tess.GetOriginalIndex(edge.neighbors.second) : tess.GetOriginalIndex(edge.neighbors.first))];
		res.insert(std::pair<size_t, ComputationalCell>(ghosts[i].second == 1 ? static_cast<size_t> (edge.neighbors.first)
			: static_cast<size_t>(edge.neighbors.second), ctemp));
	}
	return res;
}

Slope PeriodicGhostGenerator::GetGhostGradient(const Tessellation& tess,
	const vector<ComputationalCell>& /*cells*/, const vector<Slope>& gradients,
	size_t ghost_index, double /*time*/, Edge const& /*edge*/, TracerStickerNames const& /*tracerstickernames*/) const
{
	return gradients[static_cast<size_t>(tess.GetOriginalIndex(static_cast<int>(ghost_index)))];
}

