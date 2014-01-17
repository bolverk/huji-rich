#include "periodic_1d.hpp"
#include "../../misc/universal_error.hpp"

Conserved Periodic1D::CalcFlux(vector<double> const& Vertices, 
			       vector<Primitive> const& Cells,
			       RiemannSolver const& rs, 
			       vector<double> const& vertex_velocity,
			       int i) const
{
  if(i==0||i==(int)Vertices.size()-1){
    if(vertex_velocity[0]!=vertex_velocity[Vertices.size()-1])
      throw UniversalError("Vertex velocity must be the same on both sides");

    return rs.Solve(Cells[Cells.size()-1],Cells[0],
		    vertex_velocity[0]);
  }
  else
    throw UniversalError("Error in Periodic1D::CalcFlux \n Applied to bulk of grid");
}
