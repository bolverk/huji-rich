#include "spatial_reconstruction.hpp"

// Primitive gradient methods
PrimitiveGradient2D::PrimitiveGradient2D(Primitive const& xi,
	Primitive const& yi):
x(xi), y(yi) {}

PrimitiveGradient2D::PrimitiveGradient2D(void):
x(), y() {}

PrimitiveGradient2D& PrimitiveGradient2D::operator+=(PrimitiveGradient2D const& pg)
{
	x += pg.x;
	y += pg.y;
	return *this;
}

Vector2D CalcCentroid(Edge const& edge)
{
	return 0.5*(edge.GetVertex(0)+
		edge.GetVertex(1));
}

SpatialReconstruction::~SpatialReconstruction(void) {}
