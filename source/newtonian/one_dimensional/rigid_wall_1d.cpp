#include "rigid_wall_1d.hpp"
#include "../../misc/universal_error.hpp"

Conserved RigidWall1D::CalcFlux
(vector<double> const& Vertices, vector<Primitive> const& Cells,
 RiemannSolver const& rs, 
 vector<double> const& vertex_velocity,
 int i) const
{
  if(i==0){
    Primitive ghost = Cells[0];
    double vv = vertex_velocity[0];
    ghost.Velocity.x = -ghost.Velocity.x;
    return rs(ghost, Cells[0], vv);
  }
  else if(i==static_cast<int>(Vertices.size())-1){
    Primitive ghost = Cells[Vertices.size()-2];
    double vv = vertex_velocity[Vertices.size()-1];
    ghost.Velocity.x = -ghost.Velocity.x;
    return rs(Cells[Vertices.size()-2], ghost,vv);
  }
  else{
    throw UniversalError("Index inside bulk of grid");
  }
}
