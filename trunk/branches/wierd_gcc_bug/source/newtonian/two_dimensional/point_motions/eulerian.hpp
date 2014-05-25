/*! \brief Eulerian point motion scheme
  \details Sets all the components velocities of all points to be zero
*/

#ifndef EULERIAN_HPP
#define EULERIAN_HPP 1

#include "../point_motion.hpp"

//! \brief Motion scheme where all point velocities are always zero
class Eulerian: public PointMotion
{
public:
  
  Vector2D CalcVelocity
  (int /*index*/, Tessellation const* /*tessellation*/,
   vector<Primitive> const& /*primitives*/,double /*time*/);
};

#endif // EULERIAN_HPP
