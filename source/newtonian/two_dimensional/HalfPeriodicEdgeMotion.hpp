#ifndef HALF_PERIODIC_EDGE_VELOCITIES_HPP
#define HALF_PERIODIC_EDGE_VELOCITIES_HPP 1

#include "edge_velocity_calculator.hpp"

//! \brief Edge velocity calculator for half periodic domains
class HalfPeriodicEdgeVelocities : public EdgeVelocityCalculator
{
public:

	HalfPeriodicEdgeVelocities(void);

	vector<Vector2D> operator()(const Tessellation& tess,const vector<Vector2D>& point_velocities) const;
};

#endif // HALF_PERIODIC_EDGE_VELOCITIES_HPP
