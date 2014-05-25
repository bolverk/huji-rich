#ifndef CUSTOMEVOLUTION_HPP
#define CUSTOMEVOLUTION_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../../tessellation/tessellation.hpp"
#include "point_motion.hpp"
#include "spatial_reconstruction.hpp"
#include "../common/equation_of_state.hpp"
#include "../common/riemann_solver.hpp"
#include "HydroBoundaryConditions.hpp"

/*!
	\author Elad Steinberg
	\brief Abstract class for cell evolution
*/
class CustomEvolution
{
public:
	/*!
		\brief Calculates the flux across an edge
		\param tessellation The tessellation
		\param cells The primitive variables
		\param dt The sim time
		\param interpolation The interpolation scheme
		\param edge The edge
		\param facevelocity The velocity of the edge
		\param rs The Reimann Solver
		\param index The index of the custom cell
		\param boundaryconditions The external hydro boundary conditions
		\param time The sim time
		\param tracers The passive tracers
		\return The flux
	*/
	virtual Conserved CalcFlux(Tessellation const* tessellation,
		vector<Primitive> const& cells,	double dt,
		SpatialReconstruction* interpolation,Edge const& edge,
		Vector2D const& facevelocity,
		RiemannSolver const& rs,int index,
		HydroBoundaryConditions const* boundaryconditions,double time,
		vector<vector<double> > const& tracers)=0; 
		/*!
		\brief Updates the primitive variables of a cell
		\param conservedintensive The intensive consrved variables
		\param eos The equation of state
		\param cells The old primitve cells
		\param index The index of the cell to update
		\return The primitive of the cell
	*/
	virtual Primitive UpdatePrimitive(vector<Conserved> const& conservedintensive,
		EquationOfState const* eos,vector<Primitive>& cells,int index)=0;
	/*!
	\brief Virtual destructor
	*/
  virtual ~CustomEvolution(void);

  virtual bool flux_indifferent(void) const;
};

#endif
