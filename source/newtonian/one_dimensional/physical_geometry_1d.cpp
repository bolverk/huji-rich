#define _USE_MATH_DEFINES
#include <cmath>
#include "physical_geometry_1d.hpp"

PhysicalGeometry1D::~PhysicalGeometry1D(void) {}

SlabSymmetry1D::SlabSymmetry1D(void) {}

double SlabSymmetry1D::calcArea(double /*radius*/) const
{
  return 1;
}

double SlabSymmetry1D::calcVolume(double radius) const
{
  return radius;
}

CylindricalSymmetry1D::CylindricalSymmetry1D(void) {}

double CylindricalSymmetry1D::calcArea(double radius) const
{
  return 2*M_PI*radius;
}

double CylindricalSymmetry1D::calcVolume(double radius) const
{
  return M_PI*pow(radius,2);
}

SphericalSymmetry1D::SphericalSymmetry1D(void) {}

double SphericalSymmetry1D::calcArea(double radius) const
{
  return 4*M_PI*pow(radius,2);
}

double SphericalSymmetry1D::calcVolume(double radius) const
{
  return (4./3.)*M_PI*pow(radius,3);
}
