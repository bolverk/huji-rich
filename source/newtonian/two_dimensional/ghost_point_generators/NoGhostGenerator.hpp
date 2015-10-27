/*! \file NoGhostGenerator.hpp
\brief Class for not creating computational cells of ghost points
\author Elad Steinberg
*/

#ifndef NOGHOSTGENERATOR_HPP
#define NOGHOSTGENERATOR_HPP 1

#include "../GhostPointGenerator.hpp"

/*! \brief Class for creating computationalcells of ghost points for rigid walls
\author Elad Steinberg
*/
class NoGhostGenerator : public GhostPointGenerator
{
public:
	boost::container::flat_map<size_t, ComputationalCell> operator() (const Tessellation& /*tess*/,
		const vector<ComputationalCell>& /*cells*/,double /*time*/) const 
	{
		return boost::container::flat_map<size_t, ComputationalCell>();
	}

	std::pair<ComputationalCell, ComputationalCell> GetGhostGradient
	(Tessellation const& /*tess*/,
	 vector<ComputationalCell> const& /*cells*/,
	 vector<std::pair<ComputationalCell, ComputationalCell> > const& /*gradients*/,
	 size_t /*ghost_index*/, double /*time*/)const
	{
		return std::pair<ComputationalCell, ComputationalCell>();
	}
};

#endif // NOGHOSTGENERATOR_HPP
