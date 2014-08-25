/*! \file CustomEvolution.hpp
  \brief Non hydrodynamical time advance
  \author Elad Steinberg
 */

#ifndef CUSTOMEVOLUTION_HPP
#define CUSTOMEVOLUTION_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../../tessellation/tessellation.hpp"
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
	virtual Conserved CalcFlux
	(Tessellation const& tessellation,
	 vector<Primitive> const& cells,
	 double dt,
	 SpatialReconstruction& interpolation,
	 Edge const& edge,
	 Vector2D const& facevelocity,
	 RiemannSolver const& rs,int index,
	 HydroBoundaryConditions const& boundaryconditions,
	 double time,
	 vector<vector<double> > const& tracers)=0;
  /*!
	\brief Updates the primitive variables of a cell
	\param conservedintensive The intensive consrved variables
	\param eos The equation of state
	\param cells The old primitve cells
	\param index The index of the cell to update
	\param tess The tessellation
	\param time The simulation time
	\param tracers The extensive scalar tracers
	\return The primitive of the cell
	*/
	virtual Primitive UpdatePrimitive
	(vector<Conserved> const& conservedintensive,
	 EquationOfState const& eos,
	 vector<Primitive>& cells,
	 int index,
	 Tessellation const& tess,
	 double time,
	 vector<vector<double> > const& tracers)=0;

	/*!
	\brief Updates the tracer variables of a cell
	\param index The index of the cell to update
	\param tracers The extensive scalar tracers
	\param tracerchange The extensive change in the tracers
	\param cells The old primitve cells
	\param tess The tessellation
	\param time The simulation time
	\return The extensive tracers of the cell
	*/
	virtual vector<double> UpdateTracer
	(int index,vector<vector<double> >
	 const& tracers,vector<vector<double> > const& tracerchange,
	 vector<Primitive> const& cells,
	 Tessellation const& tess,
	 double time)=0;

	/*! \brief Calculates the tracer flux for a custom evolution cell
	\param tess Point and edge positions
	\param cells The hydro cells
	\param tracers The tracers
	\param dm The mass flux through the edge
	\param edge Boundary edge
	\param index The cell's index
	\param dt The time step
	\param time The sim time
	\param interp Scalar interpolation method
	\param vface The velocity of the edge
	\return The flux of the tracer
	*/
	virtual vector<double> CalcTracerFlux
	(Tessellation const& tess,
	 vector<Primitive> const& cells,vector<vector<double> > const& tracers,
	 double dm,Edge const& edge,int index,double dt,double time,
	 SpatialReconstruction const& interp,Vector2D const& vface) = 0;

	/*!
	\brief Virtual destructor
	*/
	virtual ~CustomEvolution(void);

	/*!
	\brief Sets whether to consider the cell in the time step criteria
	*/
	virtual bool TimeStepRelevant(void)const;

	  /*! \brief Calculates the velocity of the point, default is no movement
    \param index Point index
    \param tessellation Positions of the points
    \param primitives Hydrodynamic variables
    \param time The simulation time
    \return Velocity of the point
   */
  virtual Vector2D CalcVelocity(int index,
				Tessellation const& tessellation,
				vector<Primitive> const& primitives,double time);
  /*!
  \brief Determines if the cell should be considered when TracerReset is called, default is true
  \return True or false
  */
  virtual bool ShouldForceTracerReset(void)const;

  virtual bool isRelevantToInterpolation(void) const;
};

//! \brief Assigns indices to pointers to custom evolution methods
class CustomEvolutionManager
{
public:

  //! \brief Class constructor
  CustomEvolutionManager(void);

  /*! \brief Add custom evolution method to the palette
    \param cep Pointer to custom evolution method
   */
  void addCustomEvolution(CustomEvolution* cep);

  /*! \brief Returns the index of a custom evolution method
    \param cpe Pointer to custom evolution method
    \return Index
   */
  size_t getIndex(CustomEvolution* cpe) const;

  /*! \brief Converts an index to custom evolution method
    \param i Index
    \return Custom evolution method
   */
  CustomEvolution* getFunction(size_t i) const;

private:
  vector<CustomEvolution*> vtable_;
};

#endif
