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
	size_t ghostindex;
	for (size_t i = 0; i < edges.size(); ++i)
	{
		ghostindex = IsBoundaryEdge(edges[i], npoints);
		if (ghostindex == 1)
			res.push_back(std::pair<size_t, size_t>(i, 1));
		if (ghostindex == 2)
			res.push_back(std::pair<size_t, size_t>(i, 2));
	}
	return res;
}
