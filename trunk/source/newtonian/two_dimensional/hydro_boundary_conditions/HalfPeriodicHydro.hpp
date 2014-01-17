/*! \file HalfPeriodicHydro.hpp
  \brief Half Periodic Half Reflective Wall Hydro Boundary Conditions. Reflective is in the y axis
  \author Elad Steinberg
 */

#ifndef HALFPERIODICHYDRO_HPP
#define HALFPERIODICHYDRO_HPP 1

#include "../HydroBoundaryConditions.hpp"
#include "../hydrodynamics_2d.hpp"
#include "../geometric_outer_boundaries/HalfPeriodicBox.hpp"
#include "PeriodicHydro.hpp"
#include "RigidWallHydro.hpp"

//! \brief Half Periodic Half Reflective Wall Hydro Boundary Conditions. Reflective is in the y axis
class HalfPeriodicHydro: public HydroBoundaryConditions
{
public:
	/*!
	\brief Class constructor
	\param obc The outer boundary conditions
	\param rs The Riemann solver
	*/
	HalfPeriodicHydro(HalfPeriodicBox const& obc,RiemannSolver &rs);
	/*!
	\brief Class destructor
	*/
	~HalfPeriodicHydro();

	Conserved CalcFlux(Tessellation const* tessellation,
		vector<Primitive> const& cells,Vector2D const& edge_velocity,
		Edge const& edge,SpatialReconstruction const* interp,double dt,
		double time) const;

	Vector2D CalcEdgeVelocity(Tessellation const* tessellation,
		vector<Vector2D> const& point_velocities,
		Edge const& edge,double time) const;

	bool IsBoundary(Edge const& edge,Tessellation const* Data)const;

	bool IsGhostCell(int i,Tessellation const* Data) const;

	Primitive GetBoundaryPrimitive(Edge const& edge,
	  Tessellation const* Data,vector<Primitive> const& cells,double time)const;

	vector<double> GetBoundaryTracers(Edge const& edge,Tessellation const* Data,
		vector<vector<double> > const& tracers,double time)const;

  vector<double> CalcTracerFlux(Tessellation const* tessellation,
	  vector<Primitive> const& cells,
	  vector<vector<double> > const& tracers,double dm,
	  Edge const& edge,int index,double dt,
	  double time,SpatialReconstruction const* interp,
	  Vector2D const& edge_velocity) const;
private:
	HalfPeriodicBox const* obc_;
	RigidWallHydro rigid;
	PeriodicHydro periodic;

  HalfPeriodicHydro(HalfPeriodicHydro const& origin);
  HalfPeriodicHydro& operator=(HalfPeriodicHydro const& origin);
};

#endif // HALFPERIODICHYDRO_HPP
