#include "HydroBoundaryConditions.hpp"
#include "../../misc/universal_error.hpp"

Vector2D HydroBoundaryConditions::Normal(Edge const& edge,
	Tessellation const& tessellation)const
{
  if(!IsBoundary(edge,tessellation)){
    return tessellation.GetMeshPoint(edge.neighbors.second)-
      tessellation.GetMeshPoint(edge.neighbors.first);
  }
  else{
    Vector2D v;
    if(!IsGhostCell(edge.neighbors.first,tessellation)){
      v = edge.vertices.first -
	tessellation.GetMeshPoint(edge.neighbors.first);
    }
    else if(!IsGhostCell(edge.neighbors.second,tessellation)){
      v = tessellation.GetMeshPoint(edge.neighbors.second) -
	edge.vertices.first;
    }
    else{
      throw UniversalError("Both neighbors are ghosts");
    }
    Vector2D p = Parallel(edge);
    return v - p*ScalarProd(v,p)/pow(abs(p),2);
  }
}

HydroBoundaryConditions::~HydroBoundaryConditions(void) {}