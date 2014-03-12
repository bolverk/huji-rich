/*! \file hdsim.hpp
  \brief One dimensional hydrodynamic simulation
  \author Almog Yalinewich
 */

#ifndef HDSIM_1D_HPP
#define HDSIM_1D_HPP 1

#include "../common/equation_of_state.hpp"
#include "spatial_distribution1d.hpp"
#include "../common/riemann_solver.hpp"
#include "vertex_motion.hpp"
#include "spatial_reconstruction1d.hpp"
#include "boundary_conditions_1d.hpp"
#include "source_term_1d.hpp"

//! \brief Container for all hydrodynamic data
class HydroSnapshot1D
{
public:

  /*! \brief Class constructor
    \param redges Location of grid edges
    \param rcells Values of the primitive variables in hydrodynamic cells
    \param rintensive Values of the intensive conserved variables
    \param rextensive Values of the extensive conserved variables
   */
  HydroSnapshot1D(vector<double> const& redges,
		  vector<Primitive> const& rcells,
		  vector<Conserved> const& rintensive,
		  vector<Conserved> const& rextensive);

  //! \brief Position of interfaces between hydrodynamic cells
  vector<double> edges;

  //! \brief Primitive variables inside each hydrodynamic cell
  vector<Primitive> cells;

  //! \brief Intensive conserved variables of each cell
  vector<Conserved> intensive;

  //! \brief Extensive conserved variables of each cell
  vector<Conserved> extensive;
};

//! \brief Newtonian hydrodynamic simulation
class hdsim1D
{
private:

  vector<double> _Vertices;

  EquationOfState const& _eos;

  vector<Primitive> _Cells;

  vector<Conserved> _Fluxes;

  vector<double> _VertexVelocity;

  vector<Conserved> _ConservedIntensive;

  vector<Conserved> _ConservedExtensive;

  SpatialReconstruction1D const& _Interpolation;

  RiemannSolver const& _rs;

  VertexMotion const& _vm;

  BoundaryConditions1D const& _bc;

  ExternalForces1D const& force_;

  double _cfl;

  double time_;

  int cycle_;

  vector<vector<double> > tracers_intensive_;
  
  vector<vector<double> > tracers_extensive_;

public:

  /*! \brief Class constructor
    \param vertices Vertices
    \param Interpolation Interpolation method
    \param density Initial spatial density distribution
    \param pressure Initial spatial pressure distribution
    \param paravelocity Initial spatial parallel velocity distribution
    \param perpvelocity Initial spatial perpendicular velocity distribution
    \param eos Equation of state
    \param rs Riemann solver
    \param vm Vertex motion. Calculates the vertex velocities
    \param bc Boundary conditions
    \param force External force
   */
  hdsim1D(vector<double> const& vertices,
	  SpatialReconstruction1D const& Interpolation,
	  SpatialDistribution1D const& density,
	  SpatialDistribution1D const& pressure,
	  SpatialDistribution1D const& paravelocity,
	  SpatialDistribution1D const& perpvelocity,
	  EquationOfState const& eos,
	  RiemannSolver const& rs,
	  VertexMotion const& vm,
	  BoundaryConditions1D const& bc,
	  ExternalForces1D const& force);

  /*! \brief Changes the value of the Courant Friedrichs Levy coefficient
    \param cfl New value of cfl number
   */
  void overrideCFL(double cfl);

  //! \brief Advances the simulation in time
  void TimeAdvance(void);

  /*! \brief Advances the simulation in time
    \param order Order of accuracy of the scheme
   */
  void TimeAdvanceRK(int order);

  //  void AddTracer(SpatialDistribution1D const& tracer);

  // Diagnostics

  /*! \brief Returns the number of vertices
   */
  int GetVertexNo(void) const;

  /*! \brief Returns the position of a vertex
    \param i Vertex index
    \return Position of cell vertex
  */
  double GetVertexPosition(size_t i) const;

  /*! \brief Returns the number of cells
   */
  int GetCellNo(void) const;

  /*! \brief Returns the hydrodynamic variables in a cell
    \param i Cell index
    \return Primitive variables in a hydrodynamic cells
  */
  Primitive GetCell(size_t i) const;
  
  /*! \brief Returns the position of the center of the cell
  \param index Cell index
  \return Position of cell center
  */
  double GetCellCenter(size_t index) const;

  /*! \brief Returns the time of the simulation
    \return Time
   */
  double GetTime(void) const;

  /*! \brief Returns the number of times time advance was called
    \return Number of simulation time advance cycles
   */
  int GetCycle(void) const;

  /*! \brief Returns the hydrodynamic fluxes
    \return Hydrodynamic fluxes
   */
  vector<Conserved> const& getFluxes(void) const;
};

#endif // HDSIM_1D_HPP
