/*! \file eulerian.hpp
  \brief Eulerian point motion scheme
  \details Sets all the components velocities of all points to be zero
*/

#ifndef EULERIAN_HPP
#define EULERIAN_HPP 1

#include "../point_motion.hpp"

//! \brief Motion scheme where all point velocities are always zero
class Eulerian: public PointMotion
{
public:

  /*! \brief Calculates the velocity of all mesh points
    \param tess The tessellation
    \param cells Hydrodynamics cells
    \param time The simulation time
    \return Velocities of the points
  */
  vector<Vector2D> operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   double time) const;
};

#endif // EULERIAN_HPP
