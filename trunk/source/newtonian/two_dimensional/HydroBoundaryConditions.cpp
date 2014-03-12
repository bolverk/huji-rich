#include "HydroBoundaryConditions.hpp"
#include "../../misc/universal_error.hpp"

Vector2D HydroBoundaryConditions::Normal(Edge const& edge,
	Tessellation const& tessellation)const
{
  if(!IsBoundary(edge,tessellation)){
    return tessellation.GetMeshPoint(edge.GetNeighbor(1))-
      tessellation.GetMeshPoint(edge.GetNeighbor(0));
  }
  else{ 
    Vector2D v;
    if(!IsGhostCell(edge.GetNeighbor(0),tessellation)){
      v = edge.GetVertex(0) - 
	tessellation.GetMeshPoint(edge.GetNeighbor(0));
    }
    else if(!IsGhostCell(edge.GetNeighbor(1),tessellation)){
      v = tessellation.GetMeshPoint(edge.GetNeighbor(1)) - 
	edge.GetVertex(0);
    }
    else{
      throw UniversalError("Both neighbors are ghosts");
    }
    Vector2D p = Parallel(edge);
    return v - p*ScalarProd(v,p)/pow(abs(p),2);
  }
}

HydroBoundaryConditions::~HydroBoundaryConditions(void) {}
