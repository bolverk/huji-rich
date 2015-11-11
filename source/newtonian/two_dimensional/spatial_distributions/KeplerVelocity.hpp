/*! \file KeplerVelocity.hpp
  \brief Keplerian Velocity spatial distribution
  \author Elad Steinberg
*/

#ifndef KEPVEL_HPP
#define KEPVEL_HPP 1

#include "../spatial_distribution2d.hpp"
#include <cmath>

//! \brief Directions
enum Direction {xdir,ydir};

//! \brief Keplerian velocity
class KeplerVelocity : public SpatialDistribution
{
public:
	/*!
		\brief Class constructor
		\param dir The component of the velocity, can be xdir or ydir
		\param Mass The mass of the central star (G=1)
		\param center The location of the star
	*/
	KeplerVelocity(Direction dir,double Mass,Vector2D const& center=Vector2D(0,0));
	double operator()(const Vector2D& point) const;
private:
	const Direction dir_;
	const double M_;
	const Vector2D center_;
};
#endif //KEPVEL_HPP
