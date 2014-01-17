/*! \file ConstantFluxEvolution.hpp
  \brief Keeps fluxes from cell edges constant
  \author Elad Steinberg
  \todo Change name to underscore seperated lowercase words, i.e. constant_flux_evolution.hpp
 */

#ifndef COSNTANTFLUXEVOLUTION_HPP
#define COSNTANTFLUXEVOLUTION_HPP 1

#include "../CustomEvolution.hpp"
#include "../hydrodynamics_2d.hpp"

/*!
\brief Custom evolution for cells that keeps the flux constant. The momentum flux is perpindicular to the edge
\author Elad Steinberg
*/
class ConstantFluxEvolution: public CustomEvolution
{
public:
	/*!
	\brief Class constructor
	\param cell The primitive to convert to flux
	\param tracer The tracer for the tracer flux
	\param eos The equation of state
	\param entropycalc A flag stating whether the cold flows settings is on or off
	*/
	ConstantFluxEvolution(Primitive const& cell,vector<double> const& tracer,
		EquationOfState const& eos,bool entropycalc=false);
	/*!
	\brief Class destructor
	*/
	~ConstantFluxEvolution(void);

	Conserved CalcFlux(Tessellation const* tessellation,
		vector<Primitive> const& cells,	double dt,
		SpatialReconstruction* interpolation,Edge const& edge,
		Vector2D const& facevelocity,
		RiemannSolver const& rs,int index,
		HydroBoundaryConditions const* boundaryconditions,double time,
		vector<vector<double> > const& tracers); 

	Primitive UpdatePrimitive(vector<Conserved> const& conservedintensive,
		EquationOfState const* eos,vector<Primitive>& cells,int index,
		Tessellation const* tess,double time,vector<vector<double> > const& tracers);

	vector<double> UpdateTracer(int index,vector<vector<double> > const& tracers,
		vector<Primitive> const& cells,Tessellation const* tess,double time);

	vector<double> CalcTracerFlux(Tessellation const* tess,
		vector<Primitive> const& cells,vector<vector<double> > const& tracers,
		double dm,Edge const& edge,int index,double dt,double time,
		SpatialReconstruction const* interp,Vector2D const& vface);
		
private:
	const Primitive cell_;
	const vector<double> tracer_;
	const EquationOfState& eos_;
	const bool entropy_;
};

#endif //COSNTANTFLUXEVOLUTION_HPP
