#ifndef PERIODIC_EDGE_VELOCITIES_HPP
#define PERIODIC_EDGE_VELOCITIES_HPP 1

#include "edge_velocity_calculator.hpp"

//! \brief Edge velocity calculator for periodic domains
class PeriodicEdgeVelocities: public EdgeVelocityCalculator
{
public:

  PeriodicEdgeVelocities(void);

  vector<Vector2D> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities) const;
};

#endif // PERIODIC_EDGE_VELOCITIES_HPP
