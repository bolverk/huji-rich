#include "PeriodicGhostGenerator.hpp"

boost::container::flat_map<size_t, ComputationalCell> PeriodicGhostGenerator::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells) const
{
	boost::container::flat_map<size_t, ComputationalCell> res;
	vector<std::pair<size_t, size_t> > ghosts = GetOuterEdgesIndeces(tess);
	for (size_t i = 0; i < ghosts.size(); ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(ghosts[i].first));
		ComputationalCell ctemp = cells[static_cast<size_t>(ghosts[i].second == 2 ? edge.neighbors.first : edge.neighbors.second)];
		res.insert(std::pair<size_t, ComputationalCell>(ghosts[i].second == 1 ? static_cast<size_t> (edge.neighbors.first)
			: static_cast<size_t>(edge.neighbors.second), ctemp));
	}
	return res;
}

std::pair<ComputationalCell, ComputationalCell> PeriodicGhostGenerator::GetGhostGradient(const Tessellation& tess,
	const vector<ComputationalCell>& /*cells*/, const vector<std::pair<ComputationalCell, ComputationalCell> >& gradients,
	size_t ghost_index) const
{
	return gradients[static_cast<size_t>(tess.GetOriginalIndex(static_cast<int>(ghost_index)))];
}

Vector2D PeriodicGhostGenerator::GetGhostVelocity(const Tessellation& tess, const vector<ComputationalCell>& /*cells*/,
	vector<Vector2D> const& point_veolcities, size_t ghost_index, Edge const& /*edge*/)const
{
	return point_veolcities[tess.GetOriginalIndex(static_cast<int>(ghost_index))];
}
