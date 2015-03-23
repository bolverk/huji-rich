/*! \file lagrangian.hpp
  \brief Eulerian point motion scheme
  \details Sets all the components velocities of all points to be zero
*/

#ifndef LAGRANGIAN_HPP
#define LAGRANGIAN_HPP 1

#include "../point_motion.hpp"

//! \brief Motion scheme where all point velocities are equal to the material velocity
class Lagrangian: public PointMotion
{
public:

  vector<Vector2D> operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const double time) const;   
};

#endif // LAGRANGIAN_HPP
