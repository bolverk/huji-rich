#ifndef INNERRIGID_HPP
#define INNERRIGID_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../common/riemann_solver.hpp"
#include "../../tessellation/tessellation.hpp"
#include "hydrodynamics_2d.hpp"
#include "InnerBoundary.hpp"
#include <cmath>
#include "hydro_boundary_conditions/RigidWallHydro.hpp"

/*! \brief Inner Rigid Hydro Boundary Conditions
\author Elad Steinberg
*/
class InnerRigid: public InnerBoundaryConditions
{
private:
	int PointNum;
	RiemannSolver const& rs_;
	HydroBoundaryConditions *hbc_;
	RigidWallHydro rhbc;

  InnerRigid(const InnerRigid& origin);
  InnerRigid& operator=(const InnerRigid& origin);
public:
	/*!
	\brief Class constructor
	\param n The number of points in the inner boundary
	\param rs Riemann solver
	\param hbc The Outer Hydro Boundary Conditions
	*/
	InnerRigid(int n,RiemannSolver const& rs,HydroBoundaryConditions *hbc=0);
	//! \brief Class destructor
	~InnerRigid();

	Conserved CalcFlux(Tessellation const* tessellation,
		vector<Primitive> const& cells,Vector2D const& edge_velocity,
		Edge const& edge,SpatialReconstruction const* interp,double dt,
		double time) const;

	Vector2D CalcEdgeVelocity(Tessellation const* tessellation,
		vector<Vector2D> const& point_velocities,
		Edge const& edge,double time) const;

	bool IsBoundary(Edge const& edge,Tessellation const* tess)const;

	bool IsGhostCell(int i,Tessellation const* tess) const;

	int GetPointNum(void)const;

	Primitive GetBoundaryPrimitive(Edge const& edge,
		Tessellation const* tess,vector<Primitive> const& cells,
		double time)const;

	vector<double> GetBoundaryTracers(Edge const& edge,
		Tessellation const* tess,vector<vector<double> > const& tracers,
		double time)const;

	vector<double> CalcTracerFlux(Tessellation const* tessellation,
		vector<vector<double> > const& tracers,double dm,
		Edge const& edge,int index,double dt,
		double time,SpatialReconstruction const* interp)const;
	/*!
	\brief Calcualtes the flux between the real cell and the rigid counterpart
	\param tessellation The tessellation
	\param cells The primitive cells
	\param edge_velocity The velocities of the edges
	\param rs The Riemann solver
	\param edge The edge to calculate
	\param interp The spatial interpolation
	\param dt The time step
	\param ci The index of the real cell in teh edge
	\return The flux
	*/
	Conserved CalcFluxCi(Tessellation const* tessellation,
		vector<Primitive> const& cells,Vector2D const& edge_velocity,
		RiemannSolver const* rs,Edge const& edge,
		SpatialReconstruction const* interp,double dt,int ci) const;

};

#endif // INNERRIGID_HPP
