/*! \file CustomOuter.hpp
  \brief Custom mix of outer Hydro Boundary Conditions. Used with SquareBox boundary.
  \author Elad Steinberg
*/

#ifndef CUSTOMOUTER_HPP
#define CUSTOMOUTER_HPP 1

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif // _MSC_VER

#include "../HydroBoundaryConditions.hpp"

//! \brief Custom mix of outer Hydro Boundary Conditions. Used with SquareBox boundary
class CustomOuter: public HydroBoundaryConditions
{
public:
	/*!
	\brief Class constructor
	\param Left Pointer to the left HydroBoundaryConditions
	\param Right Pointer to the right HydroBoundaryConditions
	\param Down Pointer to the bottom HydroBoundaryConditions
	\param Up Pointer to the upper HydroBoundaryConditions
	*/
	CustomOuter(HydroBoundaryConditions& Left,HydroBoundaryConditions& Right,
		HydroBoundaryConditions& Down,HydroBoundaryConditions& Up);

	Conserved CalcFlux
	(Tessellation const& tessellation,
	 vector<Primitive> const& cells,
	 Vector2D const& edge_velocity,
	 Edge const& edge,SpatialReconstruction const& interp,
	 double dt, double time) const;

	Vector2D CalcEdgeVelocity(Tessellation const& tessellation,
		vector<Vector2D> const& point_velocities,
				  Edge const& edge, double time) const;

	bool IsBoundary(Edge const& edge,Tessellation const& Data)const;

	bool IsGhostCell(int i,Tessellation const& Data) const;

	Primitive GetBoundaryPrimitive(Edge const& edge,
	  Tessellation const& Data,vector<Primitive> const& cells,double
	  time)const;

	vector<double> GetBoundaryTracers(Edge const& edge,Tessellation const& Data,
		vector<vector<double> > const& tracers,double time)const;

	vector<double> CalcTracerFlux(Tessellation const& tessellation,
		vector<Primitive> const& cells,
	  vector<vector<double> > const& tracers,double dm,
	  Edge const& edge,int index,double dt,
	  double time,SpatialReconstruction const& interp,
	  Vector2D const& edge_velocity) const;
private:
  HydroBoundaryConditions& _left;
  HydroBoundaryConditions& _right;
  HydroBoundaryConditions& _down;
  HydroBoundaryConditions& _up;

  CustomOuter(const CustomOuter& origin);
  CustomOuter& operator=(const CustomOuter& origin);
};

#endif // CUSTOMOUTER_HPP
