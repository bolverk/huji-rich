#ifndef STATIONARY_BOX_HPP
#define STATIONARY_BOX_HPP 1

#include "edge_velocity_calculator.hpp"

//! \brief Edge velocity calculator for a stationary box
class StationaryBox: public EdgeVelocityCalculator
{
public:

  StationaryBox(void);

  vector<Vector2D> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities) const;
};

#endif // STATIONARY_BOX_HPP
