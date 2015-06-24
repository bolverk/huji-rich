/*! \file point_motion.hpp
  \brief Abstract class for motion of the mesh generating points
  \author Almog Yalinewich
 */

#ifndef POINT_MOTION_HPP
#define POINT_MOTION_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../../tessellation/tessellation.hpp"
#include "CustomEvolution.hpp"

/*! \brief Abstract class for motion of mesh generating points
 */
class PointMotion
{
public:

  /*! \brief Calculates the velocity of the point
    \param index Point index
    \param tessellation Positions of the points
    \param primitives Hydrodynamic variables
    \param time The simulation time
    \return Velocity of the point
   */
  virtual Vector2D CalcVelocity(int index,
				Tessellation const& tessellation,
				vector<Primitive> const& primitives,double time)= 0;

  /*! \brief Calculates the velocity of all mesh points
    \param tess The tessellation
    \param cells Hydrodynamics cells
    \param time The simulation time
	\param cevolve Custom evolution of points
	\param tracers The intensive tracers
    \return Velocities of the points
   */
  virtual vector<Vector2D> calcAllVelocities(Tessellation const& tess,
					     vector<Primitive> const& cells,
					     double time,
					     vector<CustomEvolution*> &cevolve,
					     const vector<vector<double> >& tracers);
  //! \brief Virtual destructor
  virtual ~PointMotion(void);

  /*! \brief Applies a small fix to the velocity of all mesh points once the time step is known
  \param tess The tessellation
  \param cells Hydrodynamics cells
  \param time The simulation time
  \param cevolve Custom evolution of points
  \param velocities Velocities of the points
  \param tracers The intensive tracers
  \param dt The time step
  */
  virtual void ApplyFix(Tessellation const& tess, vector<Primitive> const& cells,double time,
	  vector<CustomEvolution*> &cevolve, const vector<vector<double> >& tracers, double dt, vector<Vector2D>
	  & velocities);
};

#endif
