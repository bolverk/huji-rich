/*! \file PeriodicHydro.hpp
  \brief Periodic Wall Hydro Boundary Conditions
  \author Elad Steinberg
*/

#ifndef PERIODICHYDRO_HPP
#define PERIODICHYDRO_HPP 1

#include "../HydroBoundaryConditions.hpp"
#include "../hydrodynamics_2d.hpp"

//! \brief Periodic Wall Hydro Boundary Conditions
class PeriodicHydro: public HydroBoundaryConditions
{
public:
	//! \brief Class constructor \param rs The Riemann solver
	PeriodicHydro(RiemannSolver const& rs);
	//! \brief Class destructor
	~PeriodicHydro(void);

	Conserved CalcFlux(Tessellation const& tessellation,
		vector<Primitive> const& cells,Vector2D const& edge_velocity,
		Edge const& edge,
		SpatialReconstruction const& interp,double dt,
		double time) const;

	Vector2D CalcEdgeVelocity
		(Tessellation const& tessellation,
		vector<Vector2D> const& point_velocities,
		Edge const& edge, double time) const;

	bool IsBoundary(Edge const& edge,Tessellation const& Data)const;

	bool IsGhostCell(int i,Tessellation const& Data) const;

	Primitive GetBoundaryPrimitive
		(Edge const& edge,
		Tessellation const& Data,vector<Primitive> const& cells,double time)const;

	vector<double> GetBoundaryTracers(Edge const& edge,Tessellation const& Data,
		vector<vector<double> > const& tracers,double time)const;

	vector<double> CalcTracerFlux(Tessellation const& tessellation,
		vector<Primitive> const& cells,
		vector<vector<double> > const& tracers,double dm,
		Edge const& edge,int index,double dt,
		double time,SpatialReconstruction const& interp,
		Vector2D const& edge_velocity) const;

private:

	RiemannSolver const& rs_;
};

#endif // PERIODICHYDRO_HPP