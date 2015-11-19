#include "FreeFlowGhostGenerator.hpp"

boost::container::flat_map<size_t, ComputationalCell> FreeFlowGenerator::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells, double /*time*/) const
{
	vector<std::pair<size_t, size_t> > outer_edges=GetOuterEdgesIndeces(tess);
	boost::container::flat_map<size_t, ComputationalCell> res;
	for (size_t i = 0; i < outer_edges.size(); ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(outer_edges[i].first));
		size_t ghost_index = static_cast<size_t>(outer_edges[i].second == 1 ? edge.neighbors.first : edge.neighbors.second);
		if (tess.GetOriginalIndex(static_cast<int>(ghost_index)) < tess.GetPointNo())
		{
			int real_cell = outer_edges[i].second == 1 ? edge.neighbors.second : edge.neighbors.first;
			res[ghost_index] = cells[static_cast<size_t>(real_cell)];
			Vector2D normal = outer_edges[i].second == 1 ? tess.GetMeshPoint(edge.neighbors.first) - tess.GetMeshPoint(edge.neighbors.second)
				: tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first);
			if (ScalarProd(normal, res[ghost_index].velocity) < 0)
			{
				normal = normal / abs(normal);
				res[ghost_index].velocity -= 2 * ScalarProd(res[ghost_index].velocity, normal)*normal;
			}
		}
	}
	return res;
}

std::pair<ComputationalCell, ComputationalCell> FreeFlowGenerator::GetGhostGradient(Tessellation const& /*tess*/,
	vector<ComputationalCell> const& /*cells*/, vector<std::pair<ComputationalCell, ComputationalCell> > const& /*gradients*/,
	size_t /*ghost_index*/, double /*time*/, Edge const& /*edge*/)const
{
	return std::pair<ComputationalCell, ComputationalCell>();
}
