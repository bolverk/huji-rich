#include "RigidWallGenerator.hpp"

namespace
{
	void ReverseNormalVelocity(ComputationalCell &cell, Edge const& edge,size_t index,Tessellation const& tess)
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

std::map<size_t, ComputationalCell> RigidWallGenerator::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells) const
{
	std::map<size_t, ComputationalCell> res;
	vector<std::pair<size_t,size_t> > ghosts = GetOuterEdgesIndeces(tess);
	for (size_t i = 0; i < ghosts.size(); ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(ghosts[i].first));
		ComputationalCell ctemp = cells[static_cast<size_t>(ghosts[i].second == 2 ? edge.neighbors.first : edge.neighbors.second)];
		ReverseNormalVelocity(ctemp, edge, ghosts[i].second, tess);
		res.insert(std::pair<size_t, ComputationalCell>(ghosts[i].second == 1 ?static_cast<size_t> (edge.neighbors.first) 
			: static_cast<size_t>(edge.neighbors.second), ctemp));
	}
	return res;
}
