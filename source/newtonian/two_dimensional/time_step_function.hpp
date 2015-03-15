/*! \file time_step_function.hpp
  \author Almog Yalinewich
  \brief Abstract class for time step function
*/

#ifndef TIME_STEP_FUNCTION_HPP
#define TIME_STEP_FUNCTION_HPP 1

#include "../../tessellation/tessellation.hpp"
#include "../common/hydrodynamic_variables.hpp"
#include "../../misc/utils.hpp"
#include "CustomEvolution.hpp"

//! \brief Abstract class for time step function
class TimeStepFunction
{
public:

  /*! \brief Calculates the time step
    \param tess Tessellation
    \param cells Primitive variables
    \param point_velocities Velocities of the mesh generating points
    \param hbc Hydrodynamic boundary conditions
    \param time Time
    \param custom_evolution Lazy list of custom evolutions
    \return Time step
   */
  virtual double operator()
  (const Tessellation& tess,
   const vector<Primitive>& cells,
   const vector<Vector2D>& point_velocities,
   const HydroBoundaryConditions& hbc,
   const double time,
   const vector<CustomEvolution*>& custom_evolution) = 0;

  //! \brief Destructor
  virtual ~TimeStepFunction(void);
};

#endif // TIME_STEP_FUNCTION_HPP
