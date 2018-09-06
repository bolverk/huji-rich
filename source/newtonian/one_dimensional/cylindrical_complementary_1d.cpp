#define _USE_MATH_DEFINES
#include <cmath>
#include "cylindrical_complementary_1d.hpp"

CylindricalComplementary1D::CylindricalComplementary1D(void) {}

Conserved CylindricalComplementary1D::operator()
(const vector<double>& vertices,
 const vector<Primitive>& cells,
 size_t point,
 double /*t*/,
 double /*dt*/) const
{
  const size_t i = static_cast<size_t>(point);
  const double p = cells.at(i).Pressure;
  const double r = 0.5*(vertices.at(i)+vertices.at(i+1));
  const double volume = 
    (4./3.)*M_PI*(pow(vertices.at(i+1),3)-
		  pow(vertices.at(i),3));
  return Conserved(0,2*Vector2D(volume*(p/r),0),0);
}
