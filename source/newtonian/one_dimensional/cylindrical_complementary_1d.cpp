#ifndef __clang__
#define _USE_MATH_DEFINES
#endif // __clang__
#include <cmath>
#include "cylindrical_complementary_1d.hpp"

CylindricalComplementary1D::CylindricalComplementary1D(void) {}

Extensive CylindricalComplementary1D::operator()
  (const SimulationState1D& state,
   size_t point,
   const vector<Extensive>& /*fluxes*/,
   const PhysicalGeometry1D& /*pg*/,
   double /*t*/,
   double /*dt*/) const
{
  const vector<ComputationalCell>& cells = state.getCells();
  const vector<double>& vertices = state.getVertices();
  const size_t i = point;
  const double p = cells.at(i).pressure;
  const double r = 0.5*(vertices.at(i)+vertices.at(i+1));
  const double volume = 
    (4./3.)*M_PI*(pow(vertices.at(i+1),3)-
		  pow(vertices.at(i),3));

  Extensive res;
  res.mass = 0;
  res.momentum = 2*volume*(p/r)*Vector2D(1,0);
  res.energy = 0;
  res.tracers = vector<double>(cells.at(point).tracers.size(),0);
  return res;
}
