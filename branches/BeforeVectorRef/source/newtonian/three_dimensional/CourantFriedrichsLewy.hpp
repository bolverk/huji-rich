/*! \file CourantFriedrichsLewy.hpp
  \brief Calculates the time step according to the CFL criterion
  \author Almog Yalinewich
 */

#ifndef COURANT_FRIEDRICHS_LEWY_HPP
#define COURANT_FRIEDRICHS_LEWY_HPP 1

#include "time_step_calculator.hpp"

//! \brief Calculates the time step according to the CFL criterion
class CourantFriedrichsLewy: public TimeStepCalculator
{
public:

  /*! \brief Class constructor
    \param cfl CFL number
   */
  CourantFriedrichsLewy(double cfl);

  double operator()(const Tessellation3D& tess,
		    const vector<ComputationalCell>& cells,
		    const EquationOfState& eos) const;

private:
  const double cfl_;
};

#endif // COURANT_FRIEDRICHS_LEWY_HPP
