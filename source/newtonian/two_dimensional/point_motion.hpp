/*! \file point_motion.hpp
  \brief Abstract class for motion of the mesh generating points
  \author Almog Yalinewich
 */

#ifndef POINT_MOTION_HPP
#define POINT_MOTION_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../../tessellation/tessellation.hpp"
#include "computational_cell_2d.hpp"

/*! \brief Abstract class for motion of mesh generating points
 */
class PointMotion
{
public:

  /*! \brief Calculates the velocity of all mesh points
    \param tess The tessellation
    \param cells Hydrodynamics cells
    \param time The simulation time
	\param tracerstickernames The names of the tracers and stickers
    \return Velocities of the points
   */
  virtual vector<Vector2D> operator()(const Tessellation& tess,const vector<ComputationalCell>& cells,
	  double time,TracerStickerNames const& tracerstickernames) const = 0;

  /*! \brief Applies a small fix to the velocity of all mesh points once the time step is known
  \param tess The tessellation
  \param cells Hydrodynamics cells
  \param time The simulation time
  \param velocities Velocities of the points
  \param dt The time step
  \param tracerstickernames The names of the tracers and stickers
  \return The new velocities
  */
  virtual vector<Vector2D> ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
	  double dt, vector<Vector2D> const& velocities,TracerStickerNames const& tracerstickernames)const;

  //! \brief Virtual destructor
  virtual ~PointMotion(void);
};

#endif
