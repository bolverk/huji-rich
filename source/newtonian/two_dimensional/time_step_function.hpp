/*! \file time_step_function.hpp
  \author Almog Yalinewich
  \brief Abstract class for time step function
*/

#ifndef TIME_STEP_FUNCTION_HPP
#define TIME_STEP_FUNCTION_HPP 1

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
  double operator()(const Tessellation& tess,
		    const vector<Primitive>& cells,
		    const vector<Vector2D>& point_velocities,
		    const HydroBoundaryConditions& hbc,
		    const double time,
		    const Index2Member<CustomEvolution*>& custom_evolution);
};

#endif // TIME_STEP_FUNCTION_HPP
