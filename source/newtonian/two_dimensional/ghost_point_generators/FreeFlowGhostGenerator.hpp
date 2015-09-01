/*! \file FreeFlowGhostenerator.hpp
\brief Class for creating computationalcells of ghost points for free flow
\author Elad Steinberg
*/

#ifndef FREEFLOW_GENERATOR_HPP
#define FREEFLOW_GENERATOR_HPP 1

#include "../GhostPointGenerator.hpp"

/*! \brief Class for creating computational cells of ghost points for free flow
\author Elad Steinberg
*/
class FreeFlowGenerator : public GhostPointGenerator
{
public:
	boost::container::flat_map<size_t, ComputationalCell> operator() (const Tessellation& tess,
		const vector<ComputationalCell>& cells) const;

	std::pair<ComputationalCell, ComputationalCell> GetGhostGradient(Tessellation const& tess,
		vector<ComputationalCell> const& cells, vector<std::pair<ComputationalCell, ComputationalCell> > const& gradients,
		size_t ghost_index)const;
};

#endif // FREEFLOW_GENERATOR_HPP
