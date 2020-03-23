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
#include "cell_updater_1d.hpp"
#include "simulation_state_1d.hpp"
#include "time_step_function_1d.hpp"
#include "../two_dimensional/extensive.hpp"
#include "extensive_updater_1d.hpp"
#include "flux_calculator_1d.hpp"

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
   
   bool is_active(void) const;

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

  SimulationState1D ss_;

  EquationOfState const& eos_;

  vector<Extensive> extensives_;

  VertexMotion const& vm_;

  SourceTerm1D const& force_;

  const TimeStepFunction1D& tsf_;

  const FluxCalculator1D& fc_;

  const ExtensiveUpdater1D& eu_;
  
  const CellUpdater1D& cu_;

  double time_;

  int cycle_;

  vector<vector<double> > tracers_intensive_;
  
  vector<vector<double> > tracers_extensive_;

  ColdFlows cold_flows_;

public:

  /*! \brief Class constructor
    \param pg Physical geometry
    \param vertices Vertices
    \param density Initial spatial density distribution
    \param pressure Initial spatial pressure distribution
    \param paravelocity Initial spatial parallel velocity distribution
    \param perpvelocity Initial spatial perpendicular velocity distribution
    \param eos Equation of state
    \param vm Vertex motion. Calculates the vertex velocities
    \param force External force
  */
  /*
  hdsim1D
  (const PhysicalGeometry1D& pg,
   const vector<double>& vertices,
   const SpatialDistribution1D& density,
   const SpatialDistribution1D& pressure,
   const SpatialDistribution1D& paravelocity,
   const SpatialDistribution1D& perpvelocity,
   const EquationOfState& eos,
   const VertexMotion& vm,
   const SourceTerm1D& force,
   const TimeStepFunction1D& tsf,
   const FluxCalculator1D& fc,
   const ExtensiveUpdater1D& eu,
   const CellUpdater1D& cu);
  */

  hdsim1D
  (const PhysicalGeometry1D& pg,
   const SimulationState1D& ss,
   const EquationOfState& eos,
   const VertexMotion& vm,
   const SourceTerm1D& force,
   const TimeStepFunction1D& tsf,
   const FluxCalculator1D& fc,
   const ExtensiveUpdater1D& eu,
   const CellUpdater1D& cu);

  //! \brief Advances the simulation in time
  void TimeAdvance(void);

  //! \brief Second order time advance
  void TimeAdvance2(void);

  //  void AddTracer(SpatialDistribution1D const& tracer);

  // Scaffolding
  const vector<Primitive> getCells(void) const;

  void setCells(const vector<Primitive>& primitives);

  // Diagnostics

  /*! \brief Returns the number of vertices
    \return Number of vertices
   */
  int GetVertexNo(void) const;

  /*! \brief Returns the position of a vertex
    \param i Vertex index
    \return Position of cell vertex
  */
  double GetVertexPosition(size_t i) const;
  
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
  //  vector<Conserved> const& getFluxes(void) const;
};

#endif // HDSIM_1D_HPP
