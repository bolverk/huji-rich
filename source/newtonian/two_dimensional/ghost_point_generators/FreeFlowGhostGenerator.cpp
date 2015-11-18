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

std::pair<ComputationalCell, ComputationalCell> FreeFlowGenerator::GetGhostGradient(Tessellation const& tess,
	vector<ComputationalCell> const& cells, vector<std::pair<ComputationalCell, ComputationalCell> > const& gradients,
	size_t ghost_index, double /*time*/, Edge const& edge)const
{
	std::pair<ComputationalCell, ComputationalCell> res = gradients[static_cast<size_t>(tess.GetOriginalIndex(static_cast<int>(ghost_index)))];
	res.first.density *= -1;
	res.first.pressure *= -1;
	res.second.density *= -1;
	res.second.pressure *= -1;
	for (size_t i = 0; i < res.first.tracers.size(); ++i)
	{
		(res.first.tracers.begin() + i)->second *= -1;
		(res.second.tracers.begin() + i)->second *= -1;
	}
	Vector2D normal = edge.neighbors.first == static_cast<int>(ghost_index) ? 
		tess.GetMeshPoint(edge.neighbors.first) - tess.GetMeshPoint(edge.neighbors.second)
		: tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first);
	if (ScalarProd(normal, cells[ghost_index].velocity) > 0)
	{
		normal = normal / abs(normal);
		res.first.velocity -= 2 * ScalarProd(res.first.velocity, normal)*normal;
		res.second.velocity -= 2 * ScalarProd(res.second.velocity, normal)*normal;
	}
	return res;
}
