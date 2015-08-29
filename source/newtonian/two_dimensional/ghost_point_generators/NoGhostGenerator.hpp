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
	std::map<size_t, ComputationalCell> operator() (const Tessellation& tess,
		const vector<ComputationalCell>& cells) const 
	{
		return std::map<size_t, ComputationalCell>();
	}
};

#endif // NOGHOSTGENERATOR_HPP
