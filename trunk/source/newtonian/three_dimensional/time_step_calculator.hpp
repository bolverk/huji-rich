/*! \file time_step_calculator.hpp
  \brief Abstract class for time step calculator
  \author Almog Yalinewich
 */

#ifndef TIME_STEP_CALCULATOR_HPP
#define TIME_STEP_CALCULATOR_HPP 1

#include "../../3D/Tessellation/Tessellation3D.hpp"
#include "computational_cell.hpp"
#include "../common/equation_of_state.hpp"

class TimeStepCalculator
{
public:

  /*! \brief Calculates the time step
    \param tess Tessellation
    \param cells Computational cells
    \param eos Equation of state
   */
  virtual double operator()
  (const Tessellation3D& tess,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos) const = 0;

  //! \brief Class destructor
  virtual ~TimeStepCalculator(void);
};

#endif  // TIME_STEP_CALCULATOR_HPP
