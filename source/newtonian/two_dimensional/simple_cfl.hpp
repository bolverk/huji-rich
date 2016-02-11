/*! \file simple_cfl.hpp
  \author Almog Yalinewich
  \brief Time step based on CFL criterion
*/

#ifndef SIMPLE_CFL_HPP
#define SIMPLE_CFL_HPP 1

#include "time_step_function.hpp"

//! \brief Calculates time step according to CFL criterion
class SimpleCFL: public TimeStepFunction
{
public:

  /*! \brief Class constructor
    \param cfl CFL number
   */
  explicit SimpleCFL(const double cfl);

  double operator()(const Tessellation& tess,
		    const vector<ComputationalCell>& cells,
		    const EquationOfState& eos,
		    const vector<Vector2D>& point_velocities,
		    const double time,TracerStickerNames const& tracerstickernames) const;

private:
  const double cfl_;
};

#endif // SIMPLE_CFL_HPP
