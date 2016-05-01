/*! \file hdsim2d.hpp
  \brief Two dimensional, newtonian simulation
  \author Almog Yalinewich
*/

#ifndef HDSIM_HPP
#define HDSIM_HPP 1

#include "computational_cell_2d.hpp"
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
#include <boost/container/small_vector.hpp>
#include "computational_cell_2d.hpp"
#ifdef RICH_MPI
#include "../../mpi/ProcessorUpdate.hpp"
#endif

//! \brief Newtonian hydrodynamic simulation
class hdsim
{
private:
#ifdef RICH_MPI
	Tessellation& proctess_;
#endif

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

  TracerStickerNames tracer_sticker_names_;

  const CacheData cache_data_;

#ifdef RICH_MPI
  const ProcessorUpdate* proc_update_;
#endif

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

#ifdef RICH_MPI
   /*! \brief Returns the processor tessellation
    \return The tessellation
  */
  const Tessellation & GetProcTessellation(void)const;
#endif

   /*! \brief Class constructor
    \param tess Voronoi tessellation method
    \param proc_update Scheme for updating the locations of the processes
    \param obc Outer boundary conditions
    \param pg Physical geometry
    \param init_cond Initial conditions
    \param eos Equation of state
    \param pointmotion Motion of the mesh generating points
    \param evc Edge velocity calculator
    \param external_force External force
    \param tsf Time step function
    \param fc Flux calculator
    \param eu Extensive updater
    \param cu Cell updater
	\param tracer_sticker_names The names of the tracers and stickers
  */
#ifdef RICH_MPI
  //! \param proctess Tessellation of the processes
  //! \param proc_update Scheme for updating the position of the processes
#endif // RICH_MPI
  hdsim
	  (
#ifdef RICH_MPI
	  Tessellation& proctess,
#endif
	  Tessellation& tess,
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
   const CellUpdater& cu,
   TracerStickerNames tracer_sticker_names = TracerStickerNames()
#ifdef RICH_MPI
		  ,const ProcessorUpdate* proc_update=0
#endif
		  );

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

  //! \brief Advances the simulation in time with area fix
  void TimeAdvanceClip(void);

  //! \brief Second order time advance
  void TimeAdvance2Heun(void);

  //! \brief Second order time advance, mid point method with area fix
  void TimeAdvance2MidPointClip(void);
  //! \brief Second order time advance, mid point method
  void TimeAdvance2MidPoint(void);

  /*! \brief Change the physical geometry
    \param pg New physical geometry
   */
  void changePhysicalGeometry(const PhysicalGeometry* pg);

    /*! \brief Sets the start time
    \param t_start Start time
   */
  void setStartTime(double t_start);

  /*! \brief Sets the cycle
  \param cycle The cycle number
  */
  void setCycle(int cycle);

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

  /*!
  \brief Returns the TracerStickerNames
  \return The TracerStickerNames of the simulation
  */
  TracerStickerNames const& GetTracerStickerNames(void)const;
};

#endif
