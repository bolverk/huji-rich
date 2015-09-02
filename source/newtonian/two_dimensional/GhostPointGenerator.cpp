#include "GhostPointGenerator.hpp"

namespace
{
	size_t IsBoundaryEdge(Edge const& edge, int npoints)
	{
		if (edge.neighbors.first >= npoints)
			return 1;
		if (edge.neighbors.second >= npoints)
			return 2;
		else
			return 0;
	}
}

GhostPointGenerator::~GhostPointGenerator(void){}

vector<std::pair<size_t, size_t> > GhostPointGenerator::GetOuterEdgesIndeces(Tessellation const& tess)const
{
	vector<Edge> const& edges = tess.getAllEdges();
	vector<std::pair<size_t, size_t> > res;
	int npoints = tess.GetPointNo();
	for (size_t i = 0; i < edges.size(); ++i)
	{
	  const size_t ghostindex = IsBoundaryEdge(edges[i], npoints);
		if (ghostindex == 1)
			res.push_back(std::pair<size_t, size_t>(i, 1));
		if (ghostindex == 2)
			res.push_back(std::pair<size_t, size_t>(i, 2));
	}
	return res;
}

namespace
{
	Vector2D ReverseNormalVelocity(Vector2D const& point, Edge const& edge, size_t index, Tessellation const& tess)
	{
		Vector2D normal;
		Vector2D res(point);
		if (index == static_cast<size_t>(edge.neighbors.first))
			normal = tess.GetMeshPoint(edge.neighbors.first) - tess.GetMeshPoint(edge.neighbors.second);
		else
			normal = tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first);
		normal = normal / abs(normal);
		res -= 2 * ScalarProd(point, normal)*normal;
		return res;
	}
}

Vector2D GhostPointGenerator::GetGhostVelocity(const Tessellation& tess, const vector<ComputationalCell>& /*cells*/,
	vector<Vector2D> const& point_veolcities, size_t ghost_index,Edge const& edge)const
{

	return ReverseNormalVelocity(point_veolcities[static_cast<size_t>(tess.GetOriginalIndex(static_cast<int>(ghost_index)))],
		edge, ghost_index, tess);
		
}
