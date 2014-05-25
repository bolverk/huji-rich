/*! \brief Hydrodynamical relations 
  \author Almog Yalinewich
*/

#ifndef HYDRODYNAMICS_HPP
#define HYDRODYNAMICS_HPP 1

#include "hydrodynamic_variables.hpp"
#include "equation_of_state.hpp"

Primitive CalcPrimitive(double density, double pressure,
			Vector2D const& velocity, 
			EquationOfState const& eos);

Primitive Conserved2Primitive
(Conserved const& c, EquationOfState const& eos);

Primitive make_eos_consistent
(Primitive const& p,
 EquationOfState const& eos);

#endif // HYDRODYNAMICS_HPP
