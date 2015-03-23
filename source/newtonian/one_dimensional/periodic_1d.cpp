#include "periodic_1d.hpp"
#include "../../misc/universal_error.hpp"
#include <cmath>

namespace {
  /*
  bool effectively_zero(double x)
  {
    return std::fabs(x)<1e-14;
  }
  */
}

Conserved Periodic1D::CalcFlux(vector<double> const& Vertices, 
			       vector<Primitive> const& Cells,
			       RiemannSolver const& rs, 
			       vector<double> const& vertex_velocity,
			       int i) const
{
  if(i==0||i==static_cast<int>(Vertices.size())-1){
    /*
    if(!effectively_zero(vertex_velocity[0]-vertex_velocity[Vertices.size()-1]))
      throw UniversalError("Vertex velocity must be the same on both sides");
    */

    return rs(Cells[Cells.size()-1],Cells[0],
	      0.5*(vertex_velocity[0]+
		   vertex_velocity[Cells.size()-1]));
  }
  else
    throw UniversalError("Error in Periodic1D::CalcFlux \n Applied to bulk of grid");
}
