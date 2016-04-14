/*! \file time_step_function.hpp
  \author Almog Yalinewich
  \brief Abstract class for time step function
*/

#ifndef TIME_STEP_FUNCTION_HPP
#define TIME_STEP_FUNCTION_HPP 1

#include "../../tessellation/tessellation.hpp"
#include "../common/hydrodynamic_variables.hpp"
#include "../../misc/utils.hpp"
#include "computational_cell_2d.hpp"
#include "../common/equation_of_state.hpp"

//! \brief Abstract class for time step function
class TimeStepFunction
{
public:

  /*! \brief Calculates the time step
    \param tess Tessellation
    \param cells Primitive variables
    \param eos Equation of state
    \param point_velocities Velocities of the mesh generating points
    \param time Time
	\param tracerstickernames The names of the tracers and stickers
    \return Time step
   */
  virtual double operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const vector<Vector2D>& point_velocities,
   const double time,
	  TracerStickerNames const& tracerstickernames) const = 0;

  //! \brief Destructor
  virtual ~TimeStepFunction(void);
};

#endif // TIME_STEP_FUNCTION_HPP
