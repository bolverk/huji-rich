#ifndef COLLIDE_DENSITY_HPP
#define COLLIDE_DENSITY_HPP 1

#include "source/newtonian/two_dimensional/spatial_distribution2d.hpp"

/*!
	\brief Class that shows an example of customizing a SpatialDistribution
*/

class CollideDensity: public SpatialDistribution
{
private:
	double density1_,density2_,background_,x1_,y1_,x2_,y2_,r1_,r2_;
public:
	/*!
		\brief Class constructor
		\param density1 The density of the first blob
		\param density2 The density of the second blob
		\param background The density of the background
		\param x1 The x coordiante of the center of the first blob
		\param y1 The y coordiante of the center of the first blob
		\param x2 The x coordiante of the center of the second blob
		\param y2 The y coordiante of the center of the second blob
		\param r1 The radius of the first blob
		\param r2 The radius of the second blob
	*/
	CollideDensity(double density1,double density2,double background,
		double x1,double y1,double x2,double y2,double r1,double r2);

  double operator()(Vector2D const& point) const;
};

#endif // COLLIDE_DENSITY_HPP
