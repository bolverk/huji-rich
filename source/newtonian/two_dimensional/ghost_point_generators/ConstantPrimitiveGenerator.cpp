#include "ConstantPrimitiveGenerator.hpp"

ConstantPrimitiveGenerator::ConstantPrimitiveGenerator(ComputationalCell const& cell) : cell_(cell){}

Slope ConstantPrimitiveGenerator::GetGhostGradient(Tessellation const& tess,
	vector<ComputationalCell> const& /*cells*/, vector<Slope> const& gradients,
	size_t ghost_index, double /*time*/,Edge const& /*edge*/,TracerStickerNames const& /*tracerstickernames*/)const
{
	if (tess.GetOriginalIndex(static_cast<int>(ghost_index)) < tess.GetPointNo())
		return Slope();
	else
		return gradients[ghost_index];
}

boost::container::flat_map<size_t, ComputationalCell> ConstantPrimitiveGenerator::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells, double /*time*/,TracerStickerNames const& /*tracerstickernames*/) const
{
	vector<std::pair<size_t, size_t> > outer_edges = GetOuterEdgesIndeces(tess);
	boost::container::flat_map<size_t, ComputationalCell> res;
	for (size_t i = 0; i < outer_edges.size(); ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(outer_edges[i].first));
		size_t ghost_index = static_cast<size_t>(outer_edges[i].second == 1 ? edge.neighbors.first : edge.neighbors.second);
		if (tess.GetOriginalIndex(static_cast<int>(ghost_index)) < tess.GetPointNo())
			res[ghost_index] = cell_;
		else
			res[ghost_index] = cells[ghost_index];
	}
	return res;
}
