#ifndef BINARY_OUTFLOW_HPP
#define BINARY_OUTFLOW_HPP 1

#include "source/newtonian/two_dimensional/spatial_distribution2d.hpp"

/*!
	\brief Example of creating a custom spatial distribution for Binary outflow/inflow
*/

enum Dir{Xaxis,Yaxis};

class BinaryOutflowVelocity: public SpatialDistribution
{
private:
	double R_,val_;
	Vector2D loc_;
	Dir axis_;

public:
	/*!
	\brief Class constructor
	\param R The radius of the star
	\param val The value of the velocity
	\param loc The location of the star
	\param axis Determines if it is xvelocity or yvelocity
	*/
	BinaryOutflowVelocity(double R,double val,Vector2D loc,Dir axis);

	double operator()(Vector2D const& point) const;
};

#endif // BINARY_OUTFLOW_HPP
