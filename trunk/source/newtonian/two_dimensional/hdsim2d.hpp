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
#include "CustomEvolution.hpp"
#include "RefineStrategy.hpp"
#include "RemovalStrategy.hpp"
#include "../../misc/utils.hpp"
#include "ResetDump.hpp"
#include "../../mpi/ProcessorUpdate.hpp"


//! \brief Newtonian hydrodynamic simulation
class hdsim
{
private:

  Tessellation& _tessellation;

  Tessellation& _proctess;
  
  vector<Primitive> _cells;

  vector<Conserved> _fluxes;

  vector<Vector2D> _pointvelocity;

  vector<Vector2D> _facevelocity;

  vector<Conserved> _conservedintensive;

  vector<Conserved> _conservedextensive;

  EquationOfState const& _eos;

  RiemannSolver const& _rs;

  SpatialReconstruction& _interpolation;

  PointMotion& _pointmotion;

  HydroBoundaryConditions const& _hbc;

  OuterBoundary const& _obc;

  SourceTerm& external_force_;

  double _cfl;

  double _time;

  double _endtime;

  int cycle_;

  vector<vector<double> > tracer_;

  bool tracer_flag_,coldflows_flag_,densityfloor_;

  double as_,bs_;

  double densityMin_,pressureMin_;

  bool EntropyReCalc_;

  double _dt_external;

  ProcessorUpdate *procupdate_;

  // The purpose of these declarations is to disable copying
  hdsim(const hdsim& origin);
  
  hdsim& operator=(const hdsim& origin);

public:

	void SetProcessorMovement(ProcessorUpdate *procupdate);

  /*! \brief Returns the tessellation
    \return The tessellation
  */
  Tessellation const& GetTessellation(void)const;

   /*! \brief Returns the processor tessellation
    \return The tessellation
  */
  Tessellation const& GetProcTessellation(void)const;

  /*! 
    \brief Sets the simulation data from an external source, rebuilds the grid and resizes the tracers if needed
    \param cells The primitive cells
    \param points The mesh points
    \param time The sim time
    \param tracers The intensive tracers
  */
  void SetData(vector<Primitive> const& cells,vector<Vector2D> const& points,
	       double time,vector<vector<double> > const& tracers);
  /*!
    \brief Rearranges the simulation data according to Hibler ordering
    \param innernum The number of first points in the data not to move
  */
  void HilbertArrange(int innernum=0);

  CustomEvolutionManager custom_evolution_manager;
  vector<size_t> custom_evolution_indices;

  /*! \brief Class constructor
    \param points Initial position of the mesh generating points
    \param tessellation Voronoi tessellation method
    \param interpolation Interpolation method
    \param density Initial spatial density distribution
    \param pressure Initial spatial pressure distribution
    \param xvelocity Initial spatial distribution of the x component of the velocity
    \param yvelocity Initial spatial distribution of the y component of the velocity
    \param eos Equation of state
    \param rs Riemann solver
    \param pointmotion Motion of the mesh generating points
    \param external_force External force
    \param hbc Hydro boundary conditions
    \param obc Outer boundary conditions
    \param EntropyCalc A flag whether to recalculate the entropy during the half time step
    \param CMvalue A flag whether to give the cell value according to CM or the location of the mesh point
  */
  hdsim(vector<Vector2D> const& points,
	Tessellation& tessellation,
	SpatialReconstruction& interpolation,
	SpatialDistribution const& density,
	SpatialDistribution const& pressure,
	SpatialDistribution const& xvelocity,
	SpatialDistribution const& yvelocity,
	EquationOfState const& eos,
	RiemannSolver const& rs,
	PointMotion& pointmotion,
	SourceTerm& external_force,
	OuterBoundary const& obc,
	HydroBoundaryConditions const& hbc,
	bool EntropyCalc=false,bool CMvalue=true);
  
   /*! \brief Class constructor
    \param points Initial position of the mesh generating points
    \param tessellation Voronoi tessellation method
	\param proctess The tessellation of the proccessors
    \param interpolation Interpolation method
    \param density Initial spatial density distribution
    \param pressure Initial spatial pressure distribution
    \param xvelocity Initial spatial distribution of the x component of the velocity
    \param yvelocity Initial spatial distribution of the y component of the velocity
    \param eos Equation of state
    \param rs Riemann solver
    \param pointmotion Motion of the mesh generating points
    \param external_force External force
    \param hbc Hydro boundary conditions
    \param obc Outer boundary conditions
    \param EntropyCalc A flag whether to recalculate the entropy during the half time step
    \param CMvalue A flag whether to give the cell value according to CM or the location of the mesh point
  */
  hdsim(vector<Vector2D> const& points,
	Tessellation& tessellation,
	Tessellation& proctess,
	SpatialReconstruction& interpolation,
	SpatialDistribution const& density,
	SpatialDistribution const& pressure,
	SpatialDistribution const& xvelocity,
	SpatialDistribution const& yvelocity,
	EquationOfState const& eos,
	RiemannSolver const& rs,
	PointMotion& pointmotion,
	SourceTerm& external_force,
	OuterBoundary const& obc,
	HydroBoundaryConditions const& hbc,bool EntropyCalc=false,bool CMvalue=true);
 

  /*! \brief Class constructor from restart file
    \param dump The ResetDump file
    \param tessellation Voronoi tessellation method
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
	SpatialReconstruction& interpolation,
	EquationOfState const& eos,RiemannSolver const& rs,
	PointMotion& pointmotion,SourceTerm& external_force,
	OuterBoundary const& obc,HydroBoundaryConditions const& hbc,
	bool EntropyCalc=false);

  void load(const ResetDump& checkpoint);

  void makeCheckpoint(ResetDump& checkpoint) const;

    /*! \brief Class constructor from restart file for MPI
    \param dump The ResetDump file
    \param tessellation Voronoi tessellation method
    \param interpolation Interpolation method
    \param eos Equation of state
    \param rs Riemann solver
    \param pointmotion Motion of the mesh generating points
    \param external_force External force
    \param hbc Hydro boundary conditions
    \param obc Outer boundary conditions
    \param EntropyCalc A flag whether to recalculate the entropy during the half time step
  */
  hdsim(ResetDump const& dump,Tessellation& tessellation,Tessellation &tproc,
	SpatialReconstruction& interpolation,
	EquationOfState const& eos,RiemannSolver const& rs,
	PointMotion& pointmotion,SourceTerm& external_force,
	OuterBoundary const& obc,HydroBoundaryConditions const& hbc,
	bool EntropyCalc=false);

  /*!
    \brief Class destructor
  */
  ~hdsim(void);

  /*! \brief Overrides the value of the cfl factor
    \param cfl_new New value of the cfl factor
  */
  void SetCfl(double cfl_new);

  /*! \brief Sets the time for the sim to end exactly
    \param endtime The ending time
  */
  void SetEndTime(double endtime);
  /*!
    \brief Dictates a time step from an external source
    \param dt The time step
  */
  void SetTimeStepExternal(double dt);

  /*! \brief Advances the simulation in time
    \return Void
  */
  void TimeAdvance(void);

  /*! \brief Advances the simulation in time, second order accuracy
    \return Void
  */
  void TimeAdvance2Mid(void);
  /*! \brief Adds a tracer to the simulation
    \param tp The spatial distribution of the tracer to add
  */
  void addTracer(SpatialDistribution const& tp);
  /*! \brief Sets the Cold flow fix criteria and turns it on
    \param as The threshold for kinetic cold flows
    \param bs The threshold for external force cold flows
  */
  void SetColdFlows(double as,double bs);

  // Diagnostics

  /*! \brief Returns the total number of edges
    \return Number of edges
  */
  int GetEdgeNo(void) const;

  /*! \brief Return an edge
    \param i Edge index
    \return ith edge
  */
  Edge GetEdge(int i) const;

  /*! \brief Returns the number of cells
    \return Number of cells
  */
  int GetCellNo(void) const;

  /*! \brief Returns the primitive variables
    \param i Cell index
    \return The cell
  */
  Primitive GetCell(int i) const;

  /*! \brief Return mesh generating points
    \param i Point index
    \return Mesh generating point
  */
  Vector2D GetMeshPoint(int i) const;

  /*! \brief Returns the flux
    \param i Edge index
    \return Flux
  */
  Conserved GetFlux(int i) const;

  /*! \brief Returns the time
    \return Time
  */
  double GetTime(void) const;

  /*! \brief Returns the volume of a certain cell
    \param index Cell index
    \return Cell volume
  */
  double GetCellVolume(int index) const;

  /*! \brief Returns the number cycles
    \return Number of times TimeAdvance was called
  */
  int GetCycle(void) const;

  /*! \brief Returns the tracers vector
    \return The tracers vector
  */
  vector<vector<double> > const& getTracers(void) const;

  /*! \brief Returns the tracers vector
    \return The tracers vector
  */
  vector<vector<double> >& getTracers(void);

  /*! \brief Resets the primitive variables based on a tracer threshold
    \param alpha The tracer threshold
    \param originalD The original density distribution.
    \param originalP The original pressure distribution.
    \param originalVx The original x velocity distribution.
    \param originalVy The original y velocity distribution.
    \param tracerindex The index in the tracer to consider
  */
  void TracerReset(double alpha,SpatialDistribution const& originalD,
		   SpatialDistribution const& originalP,SpatialDistribution const& originalVx,
		   SpatialDistribution const& originalVy,int tracerindex);
 
  /*!
    \brief Returns the Courant number
    \return The courant number
  */
  double GetCfl(void)const;

  /*!
    \brief Returns the coldflow flag
  */
  bool GetColdFlowFlag(void)const;

  /*!
    \brief Returns the cold flows parameters
  */
  void GetColdFlowParm(double &a,double &b)const;

  /*!
    \brief Returns the density floor flag
  */
  bool GetDensityFloorFlag(void)const;

  /*!
    \brief Returns the density floor parameters
  */
  void GetDensityFloorParm(double &density,double &pressure)const;

  /*!
    \brief Sets a minimum density
    \param density The minimum density
    \param pressure The minimum pressure that goes along woth the minimum density
  */
  void SetDensityFloor(double density,double pressure);

  /*!
    \brief Returns all the cells
  */
  vector<Primitive>& GetAllCells(void);

  /*!
    \brief Refines cells resolution
    \param refine The refine strategy class
    \param Removed The indeces of the removed cells
    \param dr The distance from the original point to place the new point
    scaled by length of the cell. Default is 1e-4.
    \return The list of cells that were refined
  */
  vector<int> RefineCells(RefineStrategy *refine,vector<int>const& Removed,
			  double dr=1e-4);

  /*!
    \brief Removes cells
    \param remove The removal strategy.
    \return The points removed
  */
  vector<int> RemoveCells(RemovalStrategy const* remove);
  /*!
    \brief Returns the flux vector
    \return The flux vector
  */
  vector<Conserved>const& GetFluxes(void)const;
  /*!
    \brief Returns the velocity of a mesh point
    \param index The index of the mesh point
    \return The velocity of the point
  */
  Vector2D GetPointVelocity(int index)const;

  /*! \brief Returns a list of all point velocities
    \return A list of all point velocities
  */
  vector<Vector2D> const& getAllPointVelocities(void) const;
};

#endif
