#ifndef CONST_PRIMITIVE_HPP
#define CONST_PRIMITIVE_HPP 1

#include "../CustomEvolution.hpp"

class ConstantPrimitive: public CustomEvolution
{
public:

  ConstantPrimitive(Primitive const& primitive);

  bool flux_indifferent(void) const;

  Conserved CalcFlux(Tessellation const* tess,
		     vector<Primitive> const& cells, 
		     double /*dt*/,
		     SpatialReconstruction* /*interpolation*/,
		     Edge const& edge,
		     Vector2D const& facevelocity,
		     RiemannSolver const& rs,int index,
		     HydroBoundaryConditions const* /*boundaryconditions*/,
		     double /*time*/,
		     vector<vector<double> > const& /*tracers*/);

  Primitive UpdatePrimitive
  (vector<Conserved> const& /*conservedintensive*/,
   EquationOfState const* /*eos*/,
   vector<Primitive>& /*cells*/,
   int /*index*/);

private:
  const Primitive primitive_;
};

#endif // CONST_PRIMITIVE_HPP
