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

//! \brief Newtonian hydrodynamic simulation
class hdsim1D
{
private:

  const PhysicalGeometry1D& pg_;

  SimulationState1D ss_;

  const EquationOfState& eos_;

  vector<Extensive> extensives_;

  const VertexMotion& vm_;

  const SourceTerm1D& force_;

  const TimeStepFunction1D& tsf_;

  const FluxCalculator1D& fc_;

  const ExtensiveUpdater1D& eu_;
  
  const CellUpdater1D& cu_;

  double time_;

  int cycle_;

  vector<vector<double> > tracers_intensive_;
  
  vector<vector<double> > tracers_extensive_;

public:

  /*! \brief Class constructor
    \param pg Physical parameter
    \param ss Computational domain and hydro cells
    \param eos Equation of state
    \param vm Vertex motion
    \param force Source term
    \param tsf Time step function
    \param fc Flux calculator
    \param eu Extensive updater
    \param cu Cell updater
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

  void recalculateExtensives(void);

  /*! \brief Access to computational domain and hydro cells
    \return Computational domain and hydro cells
   */
  const SimulationState1D& getState(void) const;

  SimulationState1D& getState(void);

  /*! \brief Access extensive variables
    \return Extensive variables
   */
  const vector<Extensive>& getExtensives(void) const;

  /*! \brief Override extensive variables
    \param new_ext New external variables
   */
  void setExtensives(const vector<Extensive>& new_ext);

  // Diagnostics

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
