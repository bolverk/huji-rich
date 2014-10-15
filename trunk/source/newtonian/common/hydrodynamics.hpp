/*! \brief Hydrodynamical relations
  \file hydrodynamics.hpp
  \author Almog Yalinewich
*/

#ifndef HYDRODYNAMICS_HPP
#define HYDRODYNAMICS_HPP 1

#include "hydrodynamic_variables.hpp"
#include "equation_of_state.hpp"
#include "../../misc/universal_error.hpp"

/*! \brief Calculates the primitive variables
  \param density Density
  \param pressure Pressure
  \param velocity Velocity
  \param eos Equation of state
  \return Primitive variables
 */
Primitive CalcPrimitive(double density, double pressure,
			Vector2D const& velocity,
			EquationOfState const& eos);

/*! \brief Calculates the primitive variables from the conserved
  \param c Conserved variables
  \param eos Equation of state
  \return Primitive variables
 */
Primitive Conserved2Primitive
(Conserved const& c, EquationOfState const& eos);

/*! \brief Takes a set of primitive variables that do not necessarily satisfy the equation of state, and uses the density and pressure and velocity to produce a set of primitive variables that does
  \param p Primitive variables
  \param eos Equation of state
  \return Primitive variables
 */
Primitive make_eos_consistent
(Primitive const& p,
 EquationOfState const& eos);

#endif // HYDRODYNAMICS_HPP
