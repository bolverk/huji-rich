/*! \file hdsim2d.hpp
  \brief Two dimensional, newtonian simulation
  \author Almog Yalinewich
*/

#ifndef HDSIM_HPP
#define HDSIM_HPP 1

#include "../common/equation_of_state.hpp"
#include "spatial_distribution2d.hpp"
#include "../common/riemann_solver.hpp"
#include "point_motion.hpp"
#include "spatial_reconstruction.hpp"
#include "../../tessellation/tessellation.hpp"
#include "../common/hydrodynamics.hpp"
#include "SourceTerm.hpp"
#include "../../tessellation/HilbertOrder.hpp"
#include "OuterBoundary.hpp"
#include "HydroBoundaryConditions.hpp"
#include "../../misc/utils.hpp"
#include "../../misc/cached_lazy_list.hpp"
#include "../../misc/lazy_list.hpp"
#include "../../misc/universal_error.hpp"
#include "extensive.hpp"
#include "RefineStrategy.hpp"
#include "RemovalStrategy.hpp"
#include "../../misc/utils.hpp"
#include "ResetDump.hpp"
#include "physical_geometry.hpp"
#include "simple_cfl.hpp"
#include "flux_calculator_2d.hpp"
#include "cell_updater_2d.hpp"
#include "extensive_updater.hpp"
#include "cache_data.hpp"
#include "edge_velocity_calculator.hpp"

//! \brief Newtonian hydrodynamic simulation
class hdsim
{
private:

  Tessellation& tess_;

  const OuterBoundary& obc_;

  const EquationOfState& eos_;

  vector<ComputationalCell> cells_;

  vector<Extensive> extensives_;

  const PointMotion& point_motion_;

  const EdgeVelocityCalculator& edge_velocity_calculator_;

  const SourceTerm& source_;

  double time_;

  int cycle_;

  const PhysicalGeometry& pg_;

  const TimeStepFunction& tsf_;

  const FluxCalculator& fc_;

  const ExtensiveUpdater& eu_;

  const CellUpdater& cu_;

  const CacheData cache_data_;

  // The purpose of these declarations is to disable copying
  hdsim(const hdsim& origin);

  hdsim& operator=(const hdsim& origin);

public:

  /*! \brief Returns a list of all computational cells
    \return List of all computational cells
   */
  const vector<ComputationalCell>& getAllCells(void) const;

  /*! \brief Access method to manually change data
    \return Reference to computational cells
   */
  vector<ComputationalCell>& getAllCells(void);

  //! \brief Recalculates the primitives from the extensive variables
  void recalculatePrimitives(void);

  /*! \brief Recalculates extensives (in case computational cells were changed manually)
   */
  void recalculateExtensives(void);

  /*! \brief Access to physical geometry
    \return Physical geometry
   */
  const PhysicalGeometry& getPhysicalGeometry(void) const;

  /*! \brief Returns the tessellation
    \return The tessellation
  */
  const Tessellation& getTessellation(void) const;

  /*! \brief Returns the tessellation
  \return The tessellation
  */
  Tessellation& getTessellation(void);

   /*! \brief Returns the processor tessellation
    \return The tessellation
  */
  Tessellation const& GetProcTessellation(void)const;

  /*!
    \brief Rearranges the simulation data according to Hibler ordering
    \param innernum The number of first points in the data not to move
  */
  void HilbertArrange(int innernum=0);

  /*! \brief Class constructor
    \param tess Voronoi tessellation method
    \param obc Outer boundary conditions
    \param pg Physical geometry
    \param init_cond Initial conditions
    \param eos Equation of state
    \param pointmotion Motion of the mesh generating points
    \param external_force External force
    \param tsf Time step function
    \param fc Flux calculator
    \param eu Extensive updater
    \param cu Cell updater
  */
  hdsim
  (Tessellation& tess,
   const OuterBoundary& obc,
   const PhysicalGeometry& pg,
   const vector<ComputationalCell>& init_cond,
   const EquationOfState& eos,
   const PointMotion& pointmotion,
   const EdgeVelocityCalculator& evc,
   const SourceTerm& external_force,
   const TimeStepFunction& tsf,
   const FluxCalculator& fc,
   const ExtensiveUpdater& eu,
   const CellUpdater& cu);

  /*! \brief Loads reset data into simulation
    \param checkpoint Reset dump
   */
  void load(const ResetDump& checkpoint);

  /*! \brief Dumps simulation data
    \param checkpoint Reset dump where data is to be written
   */
  void makeCheckpoint(ResetDump& checkpoint) const;

  /*!
    \brief Class destructor
  */
  ~hdsim(void);

  //! \brief Advances the simulation in time
  void TimeAdvance(void);

  //! \brief Second order tiem advance
  void TimeAdvance2Heun(void);

  /*! \brief Change the physical geometry
    \param pg New physical geometry
   */
  void changePhysicalGeometry(const PhysicalGeometry* pg);

  /*! \brief Adds a tracer to the simulation
    \param name Name of tracer
    \param tp The spatial distribution of the tracer to add
  */
  void addTracer(const std::string& name,
		 const SpatialDistribution& tp);

  /*! \brief Sets the start time
    \param t_start Start time
   */
  void setStartTime(double t_start);

  // Diagnostics

  /*! \brief Returns the time
    \return Time
  */
  double getTime(void) const;

  /*! \brief Returns the volume of a certain cell
    \param index Cell index
    \return Cell volume
  */
  double getCellVolume(size_t index) const;

  /*! \brief Returns the number cycles
    \return Number of times TimeAdvance was called
  */
  int getCycle(void) const;

  /*! \brief Change the time step function
    \param tsf New time step function
   */
  void setTimeStepFunction(TimeStepFunction& tsf);

  /*! \brief Access to the equation of state
    \return Equation of state
   */
  const EquationOfState& getEos(void) const;

  /*! \brief Access to the geometric outer boundary 
    \return Geometry outer boundary
   */
  const OuterBoundary& getOuterBoundary(void) const;

  /*! \brief Returns a list of extensive variables
    \return List of extensive variables
   */
  const vector<Extensive>& getAllExtensives(void) const;

  /*! \brief Returns a list of extensive variables
    \return List of extensive variables
  */
  vector<Extensive>& getAllExtensives(void);

  /*! \brief Returns reference to the cached data
    \return Cached data
   */
  const CacheData& getCacheData(void) const;
};

#endif
