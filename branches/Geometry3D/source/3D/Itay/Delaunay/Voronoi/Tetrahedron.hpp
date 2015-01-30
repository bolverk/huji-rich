/* \file Tetrahedron.hpp
 \brief A simple encapsulation of a tetrahedron 
 \author Itay Zandbank
*/

#include <vector>
#include "GeometryCommon/Vector3D.hpp"
#include <boost/optional.hpp>

class Tetrahedron
{
private:
	std::vector<Vector3D> _vertices;
	mutable boost::optional<Vector3D> _center;   // Caches the center, which doesn't really change the object
	mutable boost::optional<double> _volume;     // Another cache

	Vector3D CalculateCenter() const;
	double CalculateVolume() const;
public:
	Tetrahedron(const std::vector<Vector3D> &vertices);
	Tetrahedron(const Vector3D v1, const Vector3D v2, const Vector3D v3, const Vector3D v4);

	Vector3D center() const;
	double volume() const;

	const vector<Vector3D> vertices() const { return _vertices;  }
	const Vector3D operator[](int index) const { return _vertices[index]; }
};

std::ostream& operator<<(std::ostream &output, const Tetrahedron &tetrahedron);