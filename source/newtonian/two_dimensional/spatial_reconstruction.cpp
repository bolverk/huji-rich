#include "spatial_reconstruction.hpp"

Vector2D CalcCentroid(Edge const& edge)
{
	return 0.5*(edge.GetVertex(0)+
		edge.GetVertex(1));
}

SpatialReconstruction::~SpatialReconstruction(void) {}
