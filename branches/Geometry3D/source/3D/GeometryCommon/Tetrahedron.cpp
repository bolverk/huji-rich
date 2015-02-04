/* \file Tetrahedron.cpp 
 \brief A simple tetrahedron wrapper
 \author Itay Zandbank */

#include "Tetrahedron.hpp"
#include "../Utilities/assert.hpp"
#include "Mat44.hpp"

using namespace std;

Tetrahedron::Tetrahedron(const std::vector<Vector3D> &vertices) : _vertices(vertices)
{
	BOOST_ASSERT(vertices.size() == 4);
}

Tetrahedron::Tetrahedron(const Vector3D v1, const Vector3D v2, const Vector3D v3, const Vector3D v4) : _vertices(4)
{
	_vertices[0] = v1;
	_vertices[1] = v2;
	_vertices[2] = v3;
	_vertices[3] = v4;
}

Vector3D Tetrahedron::center() const
{
	if (!_center.is_initialized())
		_center = CalculateCenter();
	return _center.value();
}

double Tetrahedron::volume() const
{
	if (!_volume.is_initialized())
		_volume = CalculateVolume();
	return _volume.value();
}

double Tetrahedron::radius() const
{
	if (!_radius.is_initialized())
		_radius = CalculateRadius();
	return _radius.value();
}

// \brief Find the circumcenter of a tetrahedron
// \param vertices - a vector of the 4 corners
// \returns The circumcenter
// \remark Taken from here: http://mathworld.wolfram.com/Circumsphere.html
Vector3D Tetrahedron::CalculateCenter() const
{
#define v1 _vertices[0]
#define v2 _vertices[1]
#define v3 _vertices[2]
#define v4 _vertices[3]

	Mat44<double> m_a{ v1.x, v1.y, v1.z, 1,
		v2.x, v2.y, v2.z, 1,
		v3.x, v3.y, v3.z, 1,
		v4.x, v4.y, v4.z, 1 };
	double a = m_a.determinant();

	Mat44<double> m_Dx = { abs2(v1), v1.y, v1.z, 1,
		abs2(v2), v2.y, v2.z, 1,
		abs2(v3), v3.y, v3.z, 1,
		abs2(v4), v4.y, v4.z, 1 };
	double Dx = m_Dx.determinant();

	Mat44<double> m_Dy = { abs2(v1), v1.x, v1.z, 1,
		abs2(v2), v2.x, v2.z, 1,
		abs2(v3), v3.x, v3.z, 1,
		abs2(v4), v4.x, v4.z, 1 };
	double Dy = -m_Dy.determinant();

	Mat44<double> m_Dz = { abs2(v1), v1.x, v1.y, 1,
		abs2(v2), v2.x, v2.y, 1,
		abs2(v3), v3.x, v3.y, 1,
		abs2(v4), v4.x, v4.y, 1 };
	double Dz = m_Dz.determinant();

	Mat44<double> m_c = { abs2(v1), v1.x, v1.y, v1.z,
		abs2(v2), v2.x, v2.y, v2.z,
		abs2(v3), v3.x, v3.y, v3.z,
		abs2(v4), v4.x, v4.y, v4.z };
	double c = m_c.determinant();

#undef v1
#undef v2
#undef v3
#undef v4

	return Vector3D(Dx / (2 * a), Dy / (2 * a), Dz / (2 * a));

}

double Tetrahedron::CalculateVolume() const
{
	// TODO: Calculate volume
	BOOST_ASSERT(false);  // No volume yet!
	return -1;
}

double Tetrahedron::CalculateRadius() const
{
	// The radius is the distance between the center and any of the vertices.
	return abs(center() - _vertices[0]);
}

std::ostream& operator<<(std::ostream &output, const Tetrahedron &tetrahedron)
{
	output << "{ ";
	for (int i = 0; i < 4; i++)
	{
		output << tetrahedron[i] << " ";
		if (i < 3)
			output << "- ";
	}
	output << "}";

	return output;
}