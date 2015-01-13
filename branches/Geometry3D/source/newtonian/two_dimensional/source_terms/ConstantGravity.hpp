/*! \file ConstantGravity.hpp
  \brief Constant gravity acceleration
  \author Elad Steinberg
*/

#ifndef CONSTANTGRAVITY_HPP
#define CONSTANTGRAVITY_HPP 1

#include "ConservativeForce.hpp"

//! \brief Acceleration due to constant gravity
class ConstantGravity: public Acceleration
{
public:
	/*!
	\brief Class constructor
	\param force The acceleration vector
	*/
	ConstantGravity(Vector2D const& force);

  Vector2D Calculate(Tessellation const& tess,
		     vector<Primitive> const& cells,
		     int point,
		     vector<Conserved> const& fluxes,
		     vector<Vector2D> const& point_velocity,
		     HydroBoundaryConditions const& hbc,
		     vector<vector<double> > const& tracers,
		     double time,double dt);

private:
	Vector2D force_;
};

#endif // CONSTANTGRAVITY_HPP
