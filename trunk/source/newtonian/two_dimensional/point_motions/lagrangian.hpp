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
  
  Vector2D CalcVelocity
  (int index, Tessellation const* /*tessellation*/,
   vector<Primitive> const& primitives,double time);
};

#endif // LAGRANGIAN_HPP
