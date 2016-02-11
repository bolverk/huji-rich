/*! \file PeriodicGhostGenerator.hpp
\brief Class for creating computationalcells of ghost points for periodic boundaries
\author Elad Steinberg
*/

#ifndef PERIODIC_POINT_GENERATOR_HPP
#define PERIODIC_POINT_GENERATOR_HPP 1

#include "../GhostPointGenerator.hpp"

/*! \brief Class for creating computationalcells of ghost points for periodic boundaries
\author Elad Steinberg
*/
class PeriodicGhostGenerator : public GhostPointGenerator
{
public:
	boost::container::flat_map<size_t, ComputationalCell> operator() (const Tessellation& tess,
		const vector<ComputationalCell>& cells, double time, TracerStickerNames const&
		tracerstickernames) const;

	Slope GetGhostGradient(const Tessellation& tess,
		const vector<ComputationalCell>& cells, const vector<Slope>& gradients,
		size_t ghost_index, double time, const Edge& edge, TracerStickerNames const&
		tracerstickernames) const;
};

#endif // PERIODIC_POINT_GENERATOR_HPP
