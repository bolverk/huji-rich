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
	boost::container::flat_map<size_t, ComputationalCell> operator() (const Tessellation& tess, const vector<ComputationalCell>& cells,
		double time) const;

	std::pair<ComputationalCell, ComputationalCell> GetGhostGradient(const Tessellation& tess,
		const vector<ComputationalCell>& cells,const vector<std::pair<ComputationalCell, ComputationalCell> >& gradients,
		size_t ghost_index, double time, Edge const& edge) const;
};

#endif // RIGIDWALL_GENERATOR_HPP
