#ifndef EDGE_VELOCITY_CALCULATOR_HPP
#define EDGE_VELOCITY_CALCULATOR_HPP 1

#include "../../tessellation/tessellation.hpp"

//! \brief Base class for a scheme to calculate the velocity on the edges
class EdgeVelocityCalculator
{
public:

  /*! \brief Calculates the velocity of the edges
    \param tess Tessellation
    \param point_velocities Point velocities
    \return Edge velocities
   */
  virtual vector<Vector2D> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities) const = 0;

  virtual ~EdgeVelocityCalculator(void);
};

#endif // EDGE_VELOCITY_CALCULATOR_HPP
