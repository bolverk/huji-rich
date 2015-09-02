#include "FreeFlowGhostGenerator.hpp"

boost::container::flat_map<size_t, ComputationalCell> FreeFlowGenerator::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells) const
{
	vector<std::pair<size_t, size_t> > outer_edges=GetOuterEdgesIndeces(tess);
	boost::container::flat_map<size_t, ComputationalCell> res;
	for (size_t i = 0; i < outer_edges.size(); ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(outer_edges[i].first));
		int real_cell = outer_edges[i].second == 1 ? edge.neighbors.second : edge.neighbors.first;
		/*
		res.insert(std::pair<size_t, ComputationalCell>(static_cast<size_t>(real_cell),
			cells[static_cast<size_t>(real_cell)]));
		*/
		res[static_cast<size_t>(real_cell)] = cells[static_cast<size_t>(real_cell)];
	}
	return res;
}

std::pair<ComputationalCell, ComputationalCell> FreeFlowGenerator::GetGhostGradient(Tessellation const& tess,
	vector<ComputationalCell> const& /*cells*/, vector<std::pair<ComputationalCell, ComputationalCell> > const& gradients,
	size_t ghost_index)const
{
	return gradients[static_cast<size_t>(tess.GetOriginalIndex(static_cast<int>(ghost_index)))];
}
