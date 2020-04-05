/*! \file flux_calculator_1d.hpp
  \author Almog Yalinewich
  \brief Abstract class for a flux calculator
 */

#ifndef FLUX_CALCULATOR_1D_HPP
#define FLUX_CALCULATOR_1D_HPP 1

#include <vector>
#include "../two_dimensional/extensive.hpp"
#include "simulation_state_1d.hpp"
#include "../common/equation_of_state.hpp"

using std::vector;

//! \brief Base class for a flux calculator
class FluxCalculator1D
{
public:

  /*! \brief Calculates the fluxes
    \param ss Computational grid and hydro cells
    \param vertex_velocity Velocity of the vertices
    \param eos Equation of state
    \param dt Time ste
    \return Fluxes
   */
  virtual vector<Extensive> operator()
  (const SimulationState1D& ss,
   const vector<double>& vertex_velocity,
   const EquationOfState& eos,
   const double dt) const = 0;

  virtual ~FluxCalculator1D(void);
};

#endif // FLUX_CALCULATOR_1D_HPP
