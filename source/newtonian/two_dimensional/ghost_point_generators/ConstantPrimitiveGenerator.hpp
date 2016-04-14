/*! \file ConstantPrimitiveGenerator.hpp
\brief Class for creating computationalcells of ghost points for free flow
\author Elad Steinberg
*/

#ifndef CONSTANT_PRIMITIVE_GENERATOR_HPP
#define CONSTANT_PRIMITIVE_GENERATOR_HPP 1

#include "../GhostPointGenerator.hpp"

/*! \brief Class for creating computational cells of ghost points for free flow
\author Elad Steinberg
*/
class ConstantPrimitiveGenerator : public GhostPointGenerator
{
private:
	ComputationalCell cell_;
public:
	/*!
	\brief Class constructor
	\param cell The primitive variable to be set on the boundary
	*/
  explicit ConstantPrimitiveGenerator(ComputationalCell const& cell);

	boost::container::flat_map<size_t, ComputationalCell> operator() (const Tessellation& tess,
		const vector<ComputationalCell>& cells, double time,TracerStickerNames const&
		tracerstickernames) const;

	Slope GetGhostGradient(const Tessellation& tess,const vector<ComputationalCell>& cells,const vector<Slope>& gradients,
		size_t ghost_index, double time, const Edge& edge,TracerStickerNames const& tracerstickernames)const;
};

#endif // CONSTANT_PRIMITIVE__GENERATOR_HPP
