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
#include "physical_geometry_1d.hpp"

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

//! \brief Container for cold flows data
class ColdFlows
{
public:

  ColdFlows(void);

  /*! \brief Activates the cold flows correction
    \param threshold Ratio between thermal and total energy below which the correction will be applied
   */
  void activate(double threshold);

  /*! \brief Initializes the list of entropies for each cell
    \param grid Computational grid
    \param cells Hydrodynamic cells
    \param eos Equation of state
   */
  void initializeEntropies(const vector<double>& grid,
			   const vector<Primitive>& cells,
			   const EquationOfState& eos);

  /*! \brief Advances the entropies in time
    \param fluxes Fluxes
    \param extensive Extensive conserved variables
    \param dt Time step
   */
  void advanceEntropies(const vector<Conserved>& fluxes,
			const vector<Conserved>& extensive,
			double dt);

  /*! \brief Calculates the primitive variables
    \param intensive Intensive conserved variables
    \param extensive Extensive conserved variables
    \param eos Equation of state
    \return List of primitive variables
   */
  vector<Primitive> retrieveAllPrimitive
  (const vector<Conserved>& intensive,
   const vector<Conserved>& extensive,
   const EquationOfState& eos) const;

private:
  bool active_;
  double threshold_;
  vector<double> entropies_;
};

//! \brief Newtonian hydrodynamic simulation
class hdsim1D
{
private:

  const PhysicalGeometry1D& pg_;

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

  SourceTerm1D const& force_;

  double _cfl;

  double time_;

  int cycle_;

  vector<vector<double> > tracers_intensive_;
  
  vector<vector<double> > tracers_extensive_;

  ColdFlows cold_flows_;

public:

  /*! \brief Class constructor
    \param pg Physical geometry
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
  hdsim1D
  (const PhysicalGeometry1D& pg,
   const vector<double>& vertices,
   const SpatialReconstruction1D& Interpolation,
   const SpatialDistribution1D& density,
   const SpatialDistribution1D& pressure,
   const SpatialDistribution1D& paravelocity,
   const SpatialDistribution1D& perpvelocity,
   const EquationOfState& eos,
   const RiemannSolver& rs,
   const VertexMotion& vm,
   const BoundaryConditions1D& bc,
   const SourceTerm1D& force);

  /*! \brief Changes the value of the Courant Friedrichs Levy coefficient
    \param cfl New value of cfl number
  */
  void overrideCFL(double cfl);

  //! \brief Advances the simulation in time
  void TimeAdvance(void);

  //! \brief Second order time advance
  void TimeAdvance2(void);

  /*! \brief Advances the simulation in time
    \param order Order of accuracy of the scheme
  */
  void TimeAdvanceRK(int order);
  
  /*! \brief Enables cold flows correction
    \param thres Threshold for the ratio between thermal and total energy below which the correction is applied
   */
  void enableColdFlows(double thres=1e-2);

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
