#include "Tetrahedron.hpp"

Tetrahedron::Tetrahedron()  {}

Tetrahedron::Tetrahedron(Tetrahedron const & other)
{
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
	for (int i = 0; i < 4; i++)
	{
		points[i] = other.points[i];
		neighbors[i] = other.neighbors[i];
	}
}

Tetrahedron::~Tetrahedron()
{}

Tetrahedron & Tetrahedron::operator=(Tetrahedron const & other)
{
	if (&other == this)
		return *this;
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
	for (int i = 0; i < 4; ++i)
	{
		points[i] = other.points[i];
		neighbors[i] = other.neighbors[i];
	}
	return *this;
}
