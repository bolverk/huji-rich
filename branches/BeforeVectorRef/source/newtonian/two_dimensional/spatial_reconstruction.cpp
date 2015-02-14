#include "spatial_reconstruction.hpp"

Vector2D CalcCentroid(Edge const& edge)
{
  return 0.5*(edge.vertices.first+edge.vertices.second);
}

SpatialReconstruction::~SpatialReconstruction(void) {}
