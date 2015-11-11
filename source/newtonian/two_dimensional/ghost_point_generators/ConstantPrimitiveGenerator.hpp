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
		const vector<ComputationalCell>& cells, double time) const;

  /*! \brief Calculates the gradient of a ghost cell
    \return Gradient
    \param tess Tessellation
    \param cells Computational cells
    \param gradients Gradients of the computational cells
    \param ghost_index Index of ghost point
    \param time Time
   */
	std::pair<ComputationalCell, ComputationalCell> GetGhostGradient(Tessellation const& tess,
		vector<ComputationalCell> const& cells, vector<std::pair<ComputationalCell, ComputationalCell> > const& gradients,
		size_t ghost_index, double time)const;
};

#endif // CONSTANT_PRIMITIVE__GENERATOR_HPP
