/*! \file noh_hbc.hpp
\brief Hydrodynamic boundary conditions for the noh problem
\author Elad Steinberg
*/

#ifndef NOH_HBC_HPP
#define NOH_HBC_HPP 1

#include "../../../../source/newtonian/two_dimensional/HydroBoundaryConditions.hpp"
#include "../../../../source/misc/universal_error.hpp"
//! \brief Noh problem hydro boundary conditions
class NohHBC : public HydroBoundaryConditions
{
public:
	/*! \brief Class constructor
	\param center Position of the center
	\param d0 The outer density
	\param v0 The inflow velocity
	\param p0 The outer pressure
	*/
	NohHBC(Vector2D const& center,double d0,double v0,double p0);

	Conserved CalcFlux(Tessellation const& tess,
		vector<Primitive> const& cells,
		Vector2D const& edge_velocity,
		Edge const& edge,
		SpatialReconstruction const& interp,
		double dt,
		double time) const;

	Primitive GetBoundaryPrimitive(Edge const& edge,
		Tessellation const& tess,vector<Primitive> const& cells,double time)
		const;

	bool IsGhostCell(int i,Tessellation const& tess) const;

	bool IsBoundary(Edge const& edge,Tessellation const& tess)const;

	Vector2D CalcEdgeVelocity(Tessellation const& tessellation,
		vector<Vector2D> const& point_velocities,
		Edge const& edge,
		double time) const;

	vector<double> GetBoundaryTracers(Edge const& edge,
		Tessellation const& tess,vector<vector<double> > const& tracers,double time)
		const;

	vector<double> CalcTracerFlux(Tessellation const& tessellation,
		vector<Primitive> const& cells,
		vector<vector<double> > const& tracers,double dm,
		Edge const& edge,int index,double dt,
		double time,SpatialReconstruction const& interp,Vector2D const& vface) const;

private:

	const Vector2D center_;
	const double d0_;
	const double v0_;
	const double p0_;
};

#endif // NOH_HBC_HPP
