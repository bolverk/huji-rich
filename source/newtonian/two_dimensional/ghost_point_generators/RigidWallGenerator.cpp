#include "RigidWallGenerator.hpp"

namespace
{
	void ReverseNormalVelocity(ComputationalCell &cell, Edge const& edge, size_t index, Tessellation const& tess)
	{
		Vector2D normal;
		if (index == 1)
			normal = tess.GetMeshPoint(edge.neighbors.first) - tess.GetMeshPoint(edge.neighbors.second);
		else
			normal = tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first);
		normal = normal / abs(normal);
		cell.velocity -= 2 * ScalarProd(cell.velocity, normal)*normal;
	}
}

boost::container::flat_map<size_t, ComputationalCell> RigidWallGenerator::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells, double /*time*/,TracerStickerNames const& /*tracerstickernames*/) const
{
	boost::container::flat_map<size_t, ComputationalCell> res;
	vector<std::pair<size_t, size_t> > ghosts = GetOuterEdgesIndeces(tess);
	for (size_t i = 0; i < ghosts.size(); ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(ghosts[i].first));
		size_t ghost_index = ghosts[i].second == 1 ? static_cast<size_t> (edge.neighbors.first)
			: static_cast<size_t>(edge.neighbors.second);
		if (tess.GetOriginalIndex(static_cast<int>(ghost_index)) < tess.GetPointNo())
		{
			ComputationalCell ctemp = cells[static_cast<size_t>(ghosts[i].second == 2 ? edge.neighbors.first : edge.neighbors.second)];
			ReverseNormalVelocity(ctemp, edge, ghosts[i].second, tess);
			res.insert(std::pair<size_t, ComputationalCell>(ghost_index, ctemp));
		}
		else
			res.insert(std::pair<size_t, ComputationalCell>(ghost_index, cells[ghost_index]));
	}
	return res;
}

Slope RigidWallGenerator::GetGhostGradient
(Tessellation const& tess,
	vector<ComputationalCell> const& cells,
	vector<Slope> const& /*gradients*/,
	size_t ghost_index, double /*time*/, Edge const& /*edge*/, TracerStickerNames const& /*tracerstickernames*/) const
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
