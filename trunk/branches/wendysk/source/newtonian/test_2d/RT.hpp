/*! \brief Repreduces the Vy and Pressure of the Arepo RT problem
  \author Elad Steinberg
 */

#ifndef RT_HPP
#define RT_HPP 1

#include "../two_dimensional/spatial_distribution2d.hpp"
#include <cmath>

class RT_velocityY: public SpatialDistribution
{

public:

  RT_velocityY();
  ~RT_velocityY();
  double EvalAt(Vector2D const& r) const;
};

//! \brief Spatial profile of the initial pressure for Rayleigh Taylor instability
class RT_Pressure: public SpatialDistribution
{
private:
	double g_;
	double rhou_,rhod_;
public:

  RT_Pressure(double g,double rhoup,double rhodown);
  ~RT_Pressure();
  double EvalAt(Vector2D const& r) const;
};



#endif
