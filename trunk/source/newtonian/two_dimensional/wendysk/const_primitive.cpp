#include "const_primitive.hpp"
#include "../hydrodynamics_2d.hpp"

ConstantPrimitive::ConstantPrimitive(Primitive const& primitive):
  primitive_(primitive) {}

bool ConstantPrimitive::flux_indifferent(void) const
{
  return true;
}

Conserved ConstantPrimitive::CalcFlux
(Tessellation const* tess,
 vector<Primitive> const& cells, 
 double /*dt*/,
 SpatialReconstruction* /*interpolation*/,
 Edge const& edge,
 Vector2D const& facevelocity,
 RiemannSolver const& rs,int index,
 HydroBoundaryConditions const* /*boundaryconditions*/,
 double /*time*/,
 vector<vector<double> > const& /*tracers*/)
{
  const Vector2D normal = 
    tess->GetMeshPoint(edge.GetNeighbor(1))-
    tess->GetMeshPoint(edge.GetNeighbor(0));
  const Vector2D parallel = Parallel(edge);
  Primitive left, right;
  if(edge.GetNeighbor(0)==index){
    left = primitive_;
    right = cells[edge.GetNeighbor(1)];
  }
  else if(edge.GetNeighbor(1)==index){
    left = cells[edge.GetNeighbor(0)];
    right = primitive_;
  }
  else
    throw "Something bad happened in ConstantPrimitive";
  return FluxInBulk(normal,
		    parallel,
		    left,
		    right,
		    facevelocity,
		    rs);
}

Primitive ConstantPrimitive::UpdatePrimitive
(vector<Conserved> const& /*conservedintensive*/,
 EquationOfState const* /*eos*/,
 vector<Primitive>& /*cells*/,
 int /*index*/,Tessellation const* /*tess*/,double /*time*/,
 vector<vector<double> > const& /*tracers*/)
{
  return primitive_;
}
