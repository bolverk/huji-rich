/*! \file point_motion.hpp
  \brief Abstract class for motion of the mesh generating points
  \author Almog Yalinewich
 */

#ifndef POINT_MOTION_HPP
#define POINT_MOTION_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../../tessellation/tessellation.hpp"
#include "computational_cell.hpp"

/*! \brief Abstract class for motion of mesh generating points
 */
class PointMotion
{
public:

  /*! \brief Calculates the velocity of all mesh points
    \param tess The tessellation
    \param cells Hydrodynamics cells
    \param time The simulation time
    \return Velocities of the points
   */
  virtual vector<Vector2D> operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   double time) const = 0;

  //! \brief Virtual destructor
  virtual ~PointMotion(void);
};

#endif
