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
};

#endif // RIGIDWALL_GENERATOR_HPP
