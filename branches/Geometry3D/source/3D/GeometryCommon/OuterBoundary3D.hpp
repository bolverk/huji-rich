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
private:
	Vector3D _frontUpperRight, _backLowerLeft;

public:
	//! \brief Constructs a Boundray instance
	//! \param kind Kind of boundry (rectangular or periodic)
	//! \param fromtUpperRight Front Upper Right coordinate of bounding box.
	//! \oaram backLowerLeft Back Lower Left coordinate of bounding box
	OuterBoundary3D(Vector3D frontUpperRight, Vector3D backLowerLeft);
	
	//! \brief Default constructor - a 1x1x1 cube - rectangular boundary
	// OuterBoundary3D();

	const Vector3D &FrontUpperRight() const { return _frontUpperRight; }
	const Vector3D &BackLowerLeft() const { return _backLowerLeft; }

	//\brief Returns the distance from the point, taking the subcube into account
	//\param pt The point
	//\param subcube The subcube - from '---' to '+++'.
	double distance(const Vector3D &pt, const Subcube subcube) const;

	//\brief returns the vector to the subcube. The vector is the shortest
	// to that subcube.
	//\param pt The Point
	//\param subcube The subcube
	//\returns The vector. The vector's size is distance(pt, subcube)
	Vector3D vector(const Vector3D &pt, const Subcube subcube) const;

	//\brief Creates a ghost point from pt, with the subcube
	//\param pt The point
	//\param subcube The subcube
	//\returns The ghost point.
	//\remarks How the ghost point is calculated depends on the boundary type
	virtual Vector3D ghost(const Vector3D &pt, const Subcube subcube) const = 0;

private:
	/*double distance_face(const Vector3D &pt, const Subcube subcube) const;
	double distance_edge(const Vector3D &pt, const Subcube subcube) const;
	double distance_point(const Vector3D &pt, const Subcube subcube) const; */

	Vector3D vector_face(const Vector3D &pt, const Subcube subcube) const;
	Vector3D vector_edge(const Vector3D &pt, const Subcube subcube) const;
	Vector3D vector_point(const Vector3D &pt, const Subcube subcube) const;
};

class RectangularBoundary3D : public OuterBoundary3D
{
public:
	RectangularBoundary3D(Vector3D frontUpperRight, Vector3D backLowerLeft);

	//\brief Creates a ghost point from pt, reflecting it into subcube
	//\param pt The Point
	//\param subcube The Subcube
	//\returns Pt reflected through the boundary into the subcube
	virtual Vector3D ghost(const Vector3D &pt, const Subcube subcube) const;
};
#endif // OUTERBOUNDARY3D_HPP
