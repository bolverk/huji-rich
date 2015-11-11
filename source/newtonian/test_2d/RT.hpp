/*! \file RT.hpp
  \brief Repreduces the Vy and Pressure of the Arepo RT problem
  \author Elad Steinberg
  \todo Give a more descriptive name, like rayleigh_taylor.hpp
 */

#ifndef RT_HPP
#define RT_HPP 1

#include "../two_dimensional/spatial_distribution2d.hpp"
#include <cmath>

//! \brief Y component of the velocity in Rayleigh Taylor initial conditions
class RT_velocityY: public SpatialDistribution
{

public:

  RT_velocityY();
  ~RT_velocityY();
  double operator()(const Vector2D& r) const;
};

//! \brief Spatial profile of the initial pressure for Rayleigh Taylor instability
class RT_Pressure: public SpatialDistribution
{
private:
	double g_;
	double rhou_,rhod_;
public:

  /*! \brief Class constructor
    \param g Adiabatc index
    \param rhoup Density of the upper layer
    \param rhodown Density of the lower layer
   */
  RT_Pressure(double g,double rhoup,double rhodown);

  ~RT_Pressure();

  double operator()(const Vector2D& r) const;
};



#endif
