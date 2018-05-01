#ifndef TETRAHEDRON_HPP
#define TETRAHEDRON_HPP 1

#include <array>

//points are ordered such as that the fourth point is above the plane defined by points 0 1 2 in a couter clockwise fashion
// neighbors are the tetra opposite to the triangle starting with the index of the vertice
class Tetrahedron
{
public:
	Tetrahedron();

	Tetrahedron(Tetrahedron const& other);

	~Tetrahedron();

	std::array<std::size_t, 4> points, neighbors;

	Tetrahedron& operator=(Tetrahedron const& other);
};

#endif //TETRAHEDRON_HPP