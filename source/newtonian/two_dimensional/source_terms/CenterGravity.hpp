/*! \file CenterGravity.hpp
  \brief Point source gravity force 
  \author Elad Steinberg
*/

#ifndef CENTERGRAVITY_HPP
#define CENTERGRAVITY_HPP 1

#include "ConservativeForce.hpp"

//! \brief Gravitational acceleration due to a pointlike mass
class CenterGravity: public Acceleration
{
public:
	/*!
	\brief Class constructor
	\param M The mass of the point source
	\param Rmin The softenning length
	\param center The location of the point source
	*/
	CenterGravity(double M,double Rmin,Vector2D center=Vector2D());

	Vector2D Calculate
	(Tessellation const& tess,
	 vector<Primitive> const& cells,
	 int point,vector<Conserved> const& fluxes,
	 vector<Vector2D> const& point_velocity,
	 HydroBoundaryConditions const& hbc,
	 vector<vector<double> > const& tracers,
	 double time,double dt);

private:
	double M_;
	double Rmin_;
	Vector2D _center;
};

#endif // CENTERGRAVITY_HPP
