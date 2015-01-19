/*! \file extensive_generator.hpp
  \brief Generates extensive conserved variables
  \author Almog Yalinewich
 */

#ifndef EXTENSIVE_GENERATOR_HPP
#define EXTENSIVE_GENERATOR_HPP 1

#include "conserved_3d.hpp"
#include "computational_cell.hpp"
#include "../common/equation_of_state.hpp"
#include "../../misc/utils.hpp"
#include "../../3D/GeometryCommon/Tessellation3D.hpp"

/*! \brief Calculates the extensive conserved variables
  \param cell Computational cell
  \param eos Equation of state
  \return Conserved variables
 */
Conserved3D calc_intensive(const ComputationalCell& cell,
			 const EquationOfState& eos);

//! \brief Generates a list of conserved variables
class ExtensiveGenerator: public Index2Member<Conserved3D>
{
public:

  /*! \brief Class constructor
    \param cells List of cells
    \param tess Tessellation
    \param eos Equation of state
   */
  ExtensiveGenerator(const vector<ComputationalCell>& cells,
		     const Tessellation3D& tess,
		     const EquationOfState& eos);

  size_t getLength(void) const;

  Conserved3D operator()(size_t i) const;

private:
  const vector<ComputationalCell>& cells_;
  const Tessellation3D& tess_;
  const EquationOfState& eos_;
};

#endif // EXTENSIVE_GENERATOR_HPP
