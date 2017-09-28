#include "Tetrahedron.hpp"



Tetrahedron::Tetrahedron() : points(boost::array<std::size_t, 4> ()), neighbors(boost::array<std::size_t, 4>())
{}

Tetrahedron::Tetrahedron(Tetrahedron const & other) : points(other.points),neighbors(other.neighbors)
{}

Tetrahedron::~Tetrahedron()
{}

Tetrahedron & Tetrahedron::operator=(Tetrahedron const & other)
{
	if (&other == this)
		return *this;
	points = other.points;
	neighbors = other.neighbors;
	return *this;
}
