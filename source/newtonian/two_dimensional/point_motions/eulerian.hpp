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

  vector<Vector2D> operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   double time, TracerStickerNames const& tracerstickersnames) const;
};

#endif // EULERIAN_HPP
