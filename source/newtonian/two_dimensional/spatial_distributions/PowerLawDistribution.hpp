/*! \file PowerLawDistribution.hpp
  \brief Power law spatial distribution
  \author Elad Steinberg
*/

#ifndef POWERLAW_HPP
#define POWERLAW_HPP 1

#include "../spatial_distribution2d.hpp"
#include <cmath>

//! \brief Power law spatial distribution
class PowerLawDistribution : public SpatialDistribution
{
public:
	/*!
		\brief Class constructor
		\param A The prefactor. The distribution is given as A*(x-center)^b
		\param b The exponent
		\param center The location of the star
	*/
	PowerLawDistribution(double A,double b,Vector2D const& center=Vector2D(0,0));

	double operator()(const Vector2D& point) const;
private:
	const double A_;
	const double b_;
	const Vector2D center_;
};
#endif //POWERLAW_HPP
