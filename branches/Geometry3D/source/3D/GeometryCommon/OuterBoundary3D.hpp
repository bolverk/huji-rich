/*! \file OuterBoundary3D.hpp
  \brief Outer Boundary Conditions
  \author Itat Zandbank
 */

#ifndef OUTERBOUNDARY3D_HPP
#define OUTERBOUNDARY3D_HPP 1

#include "Subcube.hpp"
//#include "Tessellation3D.hpp"

//! \brief Class describing the boundry of the Voronoi Tessallation.
class OuterBoundary3D
{
public:
	//! \brief The kind of boundry - rectangular and rigid or period.
	enum Kinds { RECTANGULAR, PERIODIC };

private:
	Vector3D _frontUpperRight, _backLowerLeft;
	Kinds _kind;

public:
	//! \brief Constructs a Boundray instance
	//! \param kind Kind of boundry (rectangular or periodic)
	//! \param fromtUpperRight Front Upper Right coordinate of bounding box.
	//! \oaram backLowerLeft Back Lower Left coordinate of bounding box
	OuterBoundary3D(Kinds kind, Vector3D frontUpperRight, Vector3D backLowerLeft);
	
	//! \brief Default constructor - a 1x1x1 cube - rectangular boundry
	OuterBoundary3D();

	const Vector3D &FrontUpperRight() const { return _frontUpperRight; }
	const Vector3D &BackLowerLeft() const { return _backLowerLeft; }

	//\brief Returns the distance from the point, taking the subcube into account
	//\param pt The point
	//\param subcube The subcube - from '---' to '+++'.
	double distance(const Vector3D &pt, const Subcube subcube) const;

private:
	double distance_face(const Vector3D &pt, const Subcube subcube) const;
	double distance_edge(const Vector3D &pt, const Subcube subcube) const;
	double distance_point(const Vector3D &pt, const Subcube subcube) const;

};
#endif // OUTERBOUNDARY3D_HPP
