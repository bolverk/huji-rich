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
	\param SoftStart The maximum radius at which the softening length works.
	\param omega The angular frequency which the source rotates with on a circular orbit
	*/
	CenterGravity(double M, double Rmin, double SoftStart = 0, Vector2D center = Vector2D(),double omega = 0);

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
	double  softlength_;
	Vector2D _center;
	double omega_;
};

#endif // CENTERGRAVITY_HPP
