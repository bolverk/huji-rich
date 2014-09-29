#ifndef CYLINDRICAL_COMPLEMENTARY_HPP
#define CYLINDRICAL_COMPLEMENTARY_HPP 1

#include "../SourceTerm.hpp"
#include "../../../tessellation/geometry.hpp"
#include "../physical_geometry.hpp"

class CylindricalComplementary: public SourceTerm
{

public:

  CylindricalComplementary(const Axis& axis);

  Conserved Calculate
  (Tessellation const& tess,
   const PhysicalGeometry& pg_,
   vector<Primitive> const& cells,
   int point,vector<Conserved> const& fluxes,
   vector<Vector2D> const& point_velocity,
   HydroBoundaryConditions const& hbc,
   vector<vector<double> > const &tracer,vector<double> &dtracer,
   vector<double> const& lengthes,double t,
   double dt);

private:
  const Axis axis_;
};


#endif // CYLINDRICAL_COMPLEMENTARY_HPP
