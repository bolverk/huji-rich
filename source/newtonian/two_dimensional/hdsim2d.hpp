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
#include "../../misc/universal_error.hpp"
#include "extensive.hpp"
#include "RefineStrategy.hpp"
#include "RemovalStrategy.hpp"
#include "../../misc/utils.hpp"
#include "ResetDump.hpp"
#include "../../mpi/ProcessorUpdate.hpp"
#include "physical_geometry.hpp"
#include "simple_cfl.hpp"
#include "flux_calculator.hpp"
#include "cell_updater.hpp"

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

  const SourceTerm& external_force_;

  double time_;

  int cycle_;

  const PhysicalGeometry& pg_;

  const TimeStepFunction& tsf_;

  const FluxCalculator& fc_;

  const CellUpdater& cu_;

  // The purpose of these declarations is to disable copying
  hdsim(const hdsim& origin);

  hdsim& operator=(const hdsim& origin);

public:

  #ifdef RICH_MPI
  /*! \brief Chooses scheme for processor motion
    \param procupdate Pointer to process motion scheme
   */
  void SetProcessorMovement(ProcessorUpdate *procupdate);
  #endif

  const vector<ComputationalCell>& getAllCells(void) const;

  const vector<Extensive>& getAllConserved(void) const;

  /*! \brief Access to physical geometry
    \return Physical geometry
   */
  const PhysicalGeometry& getPhysicalGeometry(void) const;

  /*! \brief Returns the tessellation
    \return The tessellation
  */
  const Tessellation& getTessellation(void) const;

   /*! \brief Returns the processor tessellation
    \return The tessellation
  */
  Tessellation const& GetProcTessellation(void)const;

  /*!
    \brief Rearranges the simulation data according to Hibler ordering
    \param innernum The number of first points in the data not to move
  */
  void HilbertArrange(int innernum=0);

  //! \brief container for the indices of the custom evolution
  vector<size_t> custom_evolution_indices;

  /*! \brief Class constructor
    \param tess Voronoi tessellation method
    \param eos Equation of state
    \param pointmotion Motion of the mesh generating points
    \param external_force External force
    \param obc Outer boundary conditions
  */
  hdsim(Tessellation& tess,
	const OuterBoundary& obc,
	const PhysicalGeometry& pg,
	const vector<ComputationalCell>& init_cond,
	const EquationOfState& eos,
	const PointMotion& pointmotion,
	const SourceTerm& external_force,
	const TimeStepFunction& tsf,
	const FluxCalculator& fc,
	const CellUpdater& cu);

  /*! \brief Loads reset data into simulation
    \param checkpoint Reset dump
   */
  void load(const ResetDump& checkpoint);

  /*! \brief Dumps simulation data
    \param checkpoint Reset dump where data is to be written
   */
  void makeCheckpoint(ResetDump& checkpoint) const;

    /*! \brief Class constructor from restart file for MPI
    \param dump The ResetDump file
    \param tessellation Voronoi tessellation method
    */
#ifdef RICH_MPI
  //!    \param tproc Tessellation of the processes
#endif // RICH_MPI
  /*!
    \param interpolation Interpolation method
    \param eos Equation of state
    \param rs Riemann solver
    \param pointmotion Motion of the mesh generating points
    \param external_force External force
    \param hbc Hydro boundary conditions
    \param obc Outer boundary conditions
    \param EntropyCalc A flag whether to recalculate the entropy during the half time step
  */
  hdsim(ResetDump const& dump,Tessellation& tessellation,
	#ifdef RICH_MPI
	Tessellation &tproc,
	#endif
	SpatialReconstruction& interpolation,
	EquationOfState const& eos,RiemannSolver const& rs,
	PointMotion& pointmotion,SourceTerm& external_force,
	OuterBoundary const& obc,HydroBoundaryConditions const& hbc,
	bool EntropyCalc=false);

  /*!
    \brief Class destructor
  */
  ~hdsim(void);

  //! \brief Advances the simulation in time
  void TimeAdvance(void);

  /*! \brief Change the physical geometry
    \param pg New physical geometry
   */
  void changePhysicalGeometry(const PhysicalGeometry* pg);

  /*! \brief Adds a tracer to the simulation
    \param tp The spatial distribution of the tracer to add
  */
  void addTracer(const std::string& name,
		 const SpatialDistribution& tp);

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
  double getCellVolume(int index) const;

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

  const vector<Extensive>& getExtensives(void) const;

  const vector<ComputationalCell>& getCells(void) const;
};

#endif
