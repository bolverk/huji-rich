//\field Subcube.hpp
//\brief The Subcube concept
//\remark The 3D space contains one main "subcube" defined by the boundary (it's not really a cube, but it's easy to think of it as one)
//		  and 26 adjacent "subcubes" - totaling in 27 subcubes (3 in each axis)
//\author Itay Zandbank

#ifndef SUBCUBE_H
#define SUBCUBE_H

#include "Utilities/assert.hpp"
#include "GeometryCommon/Vector3D.hpp"
#include "GeometryCommon/OuterBoundary3D.hpp"

//\brief A subcube is determined by its three offsets - one per axis. An offset can be '-', ' ' or '+' for minus, center and plus respectively
class Subcube
{
private:
	char _offsets[3];

	static size_t Num(char c)
	{
		switch (c)
		{
		case MINUS:
			return 0;
		case CENTER:
			return 1;
		case PLUS:
			return 2;
		}
		BOOST_ASSERT(false); // This value is illegal
		return 9999;
	}
public:
	Subcube(const char offsets[3]);
	const char operator[](int i) { return _offsets[i]; }

	const static char MINUS = '-';
	const static char CENTER = ' ';
	const static char PLUS = '+';

	size_t Num() const
	{
		return Num(_offsets[0]) * 9 + Num(_offsets[1]) * 3 + Num(_offsets[2]);
	}

	int NonCenters() const
	{
		return (int)(_offsets[0] != CENTER) + (int)(_offsets[1] != CENTER) + (int)(_offsets[2] != CENTER);
	}
};

bool operator<(const Subcube &sc1, const Subcube &sc2);
bool operator==(const Subcube &sc1, const Subcube &sc2);

namespace std
{
	template<>
	struct hash < Subcube >
	{
		typedef Subcube argument_type;
		typedef std::size_t result_type;

		result_type operator()(const argument_type &sc)
		{
			return sc.Num();
		}
	};
}

//\brief Calculates the distance of a point from the boundary, taking the subcube into account.
//\param pt The point
//\param boundary The boundary
//\param subcube The subcube - from '---' to '+++'.
double distance(const Vector3D &pt, const OuterBoundary3D &boundary, Subcube subcube);

#endif