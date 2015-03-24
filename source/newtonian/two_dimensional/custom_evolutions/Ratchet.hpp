/*! \file Ratchet.hpp
\brief Boundary conditions where matter is only allowed to flow in one direction
\author Almog Yalinewich
*/

#ifndef RATCHET_HPP
#define RATCHET_HPP 1

#include "../CustomEvolution.hpp"
#include "../physical_geometry.hpp"

//! \brief Cell that only allows either inflow or outflow
class Ratchet: public CustomEvolution
{
public:

	//! \brief Possible ratchet flow directions
	enum DIRECTION {in, out};

	/*! \brief Class constructor
	\param dir Direction of the flow
	\param pg Physical geometry
	*/
  Ratchet(DIRECTION dir, const PhysicalGeometry& pg);

	Conserved CalcFlux(Tessellation const& tess,
		vector<Primitive> const& cells,
		double /*dt*/,
		SpatialReconstruction& /*interp*/,
		Edge const& edge,
		Vector2D const& face_velocity,
		RiemannSolver const& rs,
		int index,
		HydroBoundaryConditions const& hbc,
		double /*time*/,
		vector<vector<double> > const& /*tracers*/);

	Primitive UpdatePrimitive(vector<Conserved> const& /*intensives*/,
		EquationOfState const& /*eos*/,vector<Primitive>& oldcells,
		int index,Tessellation const& /*tess*/,double /*time*/,
		vector<vector<double> > const& /*tracers*/);

	vector<double> UpdateTracer(int index,vector<vector<double> > const& tracers,
		vector<vector<double> > const& tracerchange,vector<Primitive> const& cells,
		Tessellation const& tess,double /*time*/);

	vector<double> CalcTracerFlux(Tessellation const& tess,
		vector<Primitive> const& cells,vector<vector<double> > const& tracers,
		double dm,Edge const& edge,int index,double dt,double time,
		SpatialReconstruction const& interp,Vector2D const& vface);

	bool TimeStepRelevant(void)const;

	bool isRelevantToInterpolation(void)const;
private:
	const DIRECTION dir_;
	const Primitive cell_;
  const PhysicalGeometry& pg_;
};

#endif // RATCHET_HPP
