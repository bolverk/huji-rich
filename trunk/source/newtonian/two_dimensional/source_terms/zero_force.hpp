/*! \file zero_force.hpp
  \brief Zero external force module
\author Elad Steinberg
*/
#ifndef ZEROFORCE_HPP
#define ZEROFORCE_HPP 1

#include "../SourceTerm.hpp"
#include "ConservativeForce.hpp"

//! \brief Zero external force module
class ZeroForce: public SourceTerm
{
public:

  Conserved Calculate
  (Tessellation const& tess,
   vector<Primitive> const& cells,
   int point,
   vector<Conserved> const& fluxes,
   vector<Vector2D> const& point_velocity,
   HydroBoundaryConditions const& hbc,
   vector<vector<double> > const& tracers,
   vector<double>& dtracer,vector<double> const& lengthes,
   double t,
   double dt);
};

//! \brief Zero acceleration
class ZeroAcceleration: public Acceleration
{
public:
  Vector2D Calculate(Tessellation const& tess,
		     vector<Primitive> const& cells,
		     int point,
		     vector<Conserved> const& fluxes,
		     vector<Vector2D> const& point_velocity,
		     HydroBoundaryConditions const& hbc,
		     vector<vector<double> > const& tracers,
		     double time,double dt);
};

#endif // ZEROFORCE_HPP