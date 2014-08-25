/*! \file RigidBodyEvolve.hpp
  \brief Cells that move like a rigid body
  \author Elad Steinberg
 */

#ifndef RIGIDBODYEVOLVE_HPP
#define RIGIDBODYEVOLVE_HPP 1

#include "../CustomEvolution.hpp"
#include "../../../misc/utils.hpp"

//! \brief Custom evolution for a cell to act as a rigid body
class RigidBodyEvolve: public CustomEvolution
{
public:
	//! \brief class constructor
	RigidBodyEvolve();

	Conserved CalcFlux(Tessellation const& tessellation,
		vector<Primitive> const& cells,	double dt,
		SpatialReconstruction& interpolation,Edge const& edge,
		Vector2D const& facevelocity,
		RiemannSolver const& rs,int index,
		HydroBoundaryConditions const& boundaryconditions,double time,
		vector<vector<double> > const& tracers);

	Primitive UpdatePrimitive(vector<Conserved> const& conservedintensive,
		EquationOfState const& eos,vector<Primitive>& cells,int index,
		Tessellation const& tess,double time,vector<vector<double> > const& tracers);

	vector<double> UpdateTracer(int index,vector<vector<double> >
		const& tracers,vector<vector<double> > const& tracerchange,vector<Primitive> const& cells,Tessellation const& tess,
		double time);

	vector<double> CalcTracerFlux(Tessellation const& tess,
		vector<Primitive> const& cells,vector<vector<double> > const& tracers,
		double dm,Edge const& edge,int index,double dt,double time,
		SpatialReconstruction const& interp,Vector2D const& vface);

	bool TimeStepRelevant(void)const;

	~RigidBodyEvolve(void);
};

#endif //RIGIDBODYEVOLVE_HPP