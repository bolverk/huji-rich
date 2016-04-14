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
		const vector<ComputationalCell>& /*cells*/,double /*time*/,TracerStickerNames const&
		/*tracerstickernames*/) const
	{
		return boost::container::flat_map<size_t, ComputationalCell>();
	}

	Slope GetGhostGradient
	(const Tessellation& /*tess*/,
	 const vector<ComputationalCell>& /*cells*/,
	 const vector<Slope>& /*gradients*/,
	 size_t /*ghost_index*/, double /*time*/, const Edge& /*edge*/,TracerStickerNames const&
		/*tracerstickernames*/) const
	{
		return Slope(ComputationalCell(), ComputationalCell());
	}
};

#endif // NOGHOSTGENERATOR_HPP
