/*! \file RigidWallGenerator.hpp
\brief Class for creating computationalcells of ghost points for rigid walls
\author Elad Steinberg
*/

#ifndef RIGIDWALL_GENERATOR_HPP
#define RIGIDWALL_GENERATOR_HPP 1

#include "../GhostPointGenerator.hpp"

/*! \brief Class for creating computationalcells of ghost points for rigid walls
\author Elad Steinberg
*/
class RigidWallGenerator : public GhostPointGenerator
{
public:
	std::map<size_t, ComputationalCell> operator() (const Tessellation& tess,
		const vector<ComputationalCell>& cells) const;

	std::pair<ComputationalCell, ComputationalCell> GetGhostGradients(Tessellation const& tess,
		vector<ComputationalCell> const& cells, vector<std::pair<ComputationalCell, ComputationalCell> > const& gradients,
		size_t ghost_index)const;

};

#endif // RIGIDWALL_GENERATOR_HPP
