/*! \file time_step_function_1d.hpp
  \author Almog Yalinewich
  \brief Abstract class for a time step function
 */

#ifndef TIME_STEP_FUNCTION_1D
#define TIME_STEP_FUNCTION_1D 1

#include "simulation_state_1d.hpp"
#include "../common/equation_of_state.hpp"

//! \brief Base class for a time step calculator
class TimeStepFunction1D
{
public:

  /*! \brief Calculates the time step
    \param ss Computational domain and hydro cells
    \param eos Equation of state
    \return Time step
   */
  virtual double operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos) const = 0;

  virtual ~TimeStepFunction1D(void);
};

#endif // TIME_STEP_FUNCTION_1D
