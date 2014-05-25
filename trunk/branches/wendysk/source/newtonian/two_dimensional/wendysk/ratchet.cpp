#include "ratchet.hpp"
#include "../hydrodynamics_2d.hpp"

Ratchet::Ratchet(DIRECTION dir,
		 vector<Primitive> const& init_cells):
  dir_(dir),
  init_cells_(init_cells) {}

bool Ratchet::flux_indifferent(void) const
{
  return true;
}

Conserved Ratchet::CalcFlux
(Tessellation const* tess,
 vector<Primitive> const& cells,
 double /*dt*/,
 SpatialReconstruction* /*interp*/,
 Edge const& edge,
 Vector2D const& face_velocity,
 RiemannSolver const& rs,
 int index,
 HydroBoundaryConditions const* hbc,
 double /*time*/,
 vector<vector<double> > const& /*tracers*/)
{
  int other;
  int my_index;
  if(edge.GetNeighbor(0)==index){
    my_index = 0;
    other = edge.GetNeighbor(1);
  }
  else{
    my_index = 1;
    other = edge.GetNeighbor(0);
  }
  if(hbc->IsGhostCell(other,tess))
    return Conserved();
  const Vector2D p = Parallel(edge);
  const Vector2D n = 
    tess->GetMeshPoint(edge.GetNeighbor(1))
    -tess->GetMeshPoint(edge.GetNeighbor(0));
  const Vector2D outward = 
    tess->GetMeshPoint(edge.GetNeighbor(1-my_index))
    -tess->GetMeshPoint(edge.GetNeighbor(my_index));
  Primitive ghost = cells[other];
  if (((dir_==in)&&(ScalarProd(ghost.Velocity,outward)>0))||
      ((dir_==out)&&(ScalarProd(ghost.Velocity,outward)<0)))
    ghost.Velocity = Reflect(cells[other].Velocity, p);
  Primitive left, right;
  if(0==my_index){
    left = ghost;
    right = cells[other];
  }
  else{
    left = cells[other];
    right = ghost;
  }
  return FluxInBulk(n,p,left,right,face_velocity,rs);
}

Primitive Ratchet::UpdatePrimitive(vector<Conserved> const& /*intensives*/,
				   EquationOfState const* /*eos*/,
				   vector<Primitive>& /*cells*/,
				   int index)
{
  return init_cells_[index];
}
