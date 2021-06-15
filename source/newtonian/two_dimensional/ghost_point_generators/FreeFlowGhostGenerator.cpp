#include "FreeFlowGhostGenerator.hpp"

boost::container::flat_map<size_t, ComputationalCell> FreeFlowGenerator::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells, double /*time*/) const
{
	vector<std::pair<size_t, size_t> > outer_edges = GetOuterEdgesIndeces(tess);
	boost::container::flat_map<size_t, ComputationalCell> res;
	for (size_t i = 0; i < outer_edges.size(); ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(outer_edges[i].first));
		size_t ghost_index = static_cast<size_t>(outer_edges[i].second == 1 ? edge.neighbors.first : edge.neighbors.second);
		if (tess.GetOriginalIndex(static_cast<int>(ghost_index)) < tess.GetPointNo())
		{
			int real_cell = outer_edges[i].second == 1 ? edge.neighbors.second : edge.neighbors.first;
			res[ghost_index] = cells[static_cast<size_t>(real_cell)];
		}
		else
			res.insert(std::pair<size_t, ComputationalCell>(ghost_index, cells[ghost_index]));
	}
	return res;
}

Slope FreeFlowGenerator::GetGhostGradient(Tessellation const& tess,
	vector<ComputationalCell> const& cells, vector<Slope> const& /*gradients*/,
	size_t ghost_index, double /*time*/, Edge const& /*edge*/)const
{
	ComputationalCell cell(cells[static_cast<size_t>(tess.GetOriginalIndex(static_cast<int>(ghost_index)))]);
	cell.density = 0;
	cell.pressure = 0;
	cell.velocity = Vector2D(0, 0);
	cell.tracers = cells[0].tracers;
	size_t N = cell.tracers.size();
	for (size_t i = 0; i < N; ++i)
		cell.tracers[i] = 0;
	return Slope(cell, cell);
}
