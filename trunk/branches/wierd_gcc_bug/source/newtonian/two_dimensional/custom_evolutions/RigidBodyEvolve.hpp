#ifndef RIGIDBODYEVOLVE_HPP
#define RIGIDBODYEVOLVE_HPP 1

#include "../CustomEvolution.hpp"

//! \brief Custom evolution for a cell to act as a rigid body
class RigidBodyEvolve: public CustomEvolution
{

public:
	//! \brief class constructor
	RigidBodyEvolve();

	Conserved CalcFlux(Tessellation const* tessellation,
		vector<Primitive> const& cells,	double dt,
		SpatialReconstruction* interpolation,Edge const& edge,
		Vector2D const& facevelocity,
		RiemannSolver const& rs,int index,
		HydroBoundaryConditions const* boundaryconditions,double time,
		vector<vector<double> > const& tracers); 

	Primitive UpdatePrimitive(vector<Conserved> const& conservedintensive,
		EquationOfState const* eos,vector<Primitive>& cells,int index);
};

#endif //RIGIDBODYEVOLVE_HPP
