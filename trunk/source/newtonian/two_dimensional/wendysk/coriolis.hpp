/*! \file coriolis.hpp
  \brief Coriolis force
  \author Almog Yalinewich
 */ 

#ifndef CORIOLIS_HPP
#define CORIOLIS_HPP 1

#include "../SourceTerm.hpp"

//! \brief Source term due to Coriolis force
class Coriolis: public SourceTerm
{
public:

  /*! \brief Class constructor
    \param angular_velocity Angular velocity
   */
  Coriolis(double angular_velocity);

  Conserved Calculate(Tessellation const& tess,
		      vector<Primitive> const& cells,
		      int point,
		      vector<Conserved> const& /*fluxes*/,
		      vector<Vector2D> const& /*point_velocity*/,
		      HydroBoundaryConditions const& /*hbc*/,
		      vector<vector<double> > const& /*tracers*/,
		      vector<double>& /*dtracer*/,
		      double /*time*/,
		      double /*dt*/);

private:

    const double omega_;
};

#endif //CORIOLIS_HPP
