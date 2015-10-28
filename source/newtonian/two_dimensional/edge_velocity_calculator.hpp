#ifndef EDGE_VELOCITY_CALCULATOR_HPP
#define EDGE_VELOCITY_CALCULATOR_HPP 1

#include "../../tessellation/tessellation.hpp"

class EdgeVelocityCalculator
{
public:

  virtual vector<Vector2D> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities) const = 0;

  virtual ~EdgeVelocityCalculator(void);
};

#endif // EDGE_VELOCITY_CALCULATOR_HPP
