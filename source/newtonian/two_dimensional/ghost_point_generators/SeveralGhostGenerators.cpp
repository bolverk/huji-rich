#include "SeveralGhostGenerators.hpp"

typedef boost::container::flat_map<size_t, ComputationalCell> GhostCells;

SeveralGhostGenerators::SeveralGhostGenerators(vector<GhostPointGenerator*> ghosts,GhostCriteria const& ghostchooser) :
ghosts_(ghosts),ghost_chooser_(ghostchooser)
{}

boost::container::flat_map<size_t, ComputationalCell> SeveralGhostGenerators::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells) const
{
	size_t nghosts = ghosts_.size();
	vector<GhostCells> ghost_cells(nghosts);
	for (size_t i = 0; i < nghosts; ++i)
		ghost_cells[i] = ghosts_[i]->operator()(tess, cells);
	GhostCells res;
	for (GhostCells::const_iterator it = ghost_cells[0].begin(); it != ghost_cells[0].end();++it)
	{
		res.insert(std::pair<size_t,ComputationalCell>(it->first,
			ghost_cells[ghost_chooser_.GhostChoose(tess, it->first)][it->first]));
	}
	return res;
}

std::pair<ComputationalCell, ComputationalCell> SeveralGhostGenerators::GetGhostGradient(const Tessellation& tess,
	const vector<ComputationalCell>& cells, const vector<std::pair<ComputationalCell, ComputationalCell> >& gradients,
	size_t ghost_index) const
{
	return ghosts_[ghost_chooser_.GhostChoose(tess, ghost_index)]->GetGhostGradient(tess, cells, gradients, ghost_index);
}

Vector2D SeveralGhostGenerators::GetGhostVelocity(const Tessellation& tess, const vector<ComputationalCell>& cells,
	vector<Vector2D> const& point_veolcities, size_t ghost_index, Edge const& edge)const
{
	return ghosts_[ghost_chooser_.GhostChoose(tess, ghost_index)]->GetGhostVelocity(tess, cells, point_veolcities,
		ghost_index,edge);
}
