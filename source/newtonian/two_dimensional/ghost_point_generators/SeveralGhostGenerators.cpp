#include "SeveralGhostGenerators.hpp"

GhostCriteria::~GhostCriteria(void) {}

typedef boost::container::flat_map<size_t, ComputationalCell> GhostCells;

SeveralGhostGenerators::SeveralGhostGenerators(vector<GhostPointGenerator*> ghosts, GhostCriteria const& ghostchooser) :
	ghosts_(ghosts), ghost_chooser_(ghostchooser)
{}

boost::container::flat_map<size_t, ComputationalCell> SeveralGhostGenerators::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells, double time, TracerStickerNames const&
	tracerstickernames) const
{
	size_t nghosts = ghosts_.size();
	vector<GhostCells> ghost_cells(nghosts);
	for (size_t i = 0; i < nghosts; ++i)
		ghost_cells[i] = ghosts_[i]->operator()(tess, cells, time, tracerstickernames);
	GhostCells res;
	for (GhostCells::const_iterator it = ghost_cells[0].begin(); it != ghost_cells[0].end(); ++it)
	{
		res.insert(std::pair<size_t, ComputationalCell>(it->first,
			ghost_cells[ghost_chooser_.GhostChoose(tess, static_cast<int>(it->first))][it->first]));
	}
	return res;
}

Slope SeveralGhostGenerators::GetGhostGradient(const Tessellation& tess,
	const vector<ComputationalCell>& cells, const vector<Slope>& gradients,
	size_t ghost_index, double time, Edge const& edge,TracerStickerNames const&
	tracerstickernames) const
{
	return ghosts_[ghost_chooser_.GhostChoose(tess, 
		static_cast<int>(ghost_index))]->GetGhostGradient(tess, cells, gradients, ghost_index, time, edge,tracerstickernames);
}

