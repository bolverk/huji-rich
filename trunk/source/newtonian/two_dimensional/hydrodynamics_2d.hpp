/*! \file hydrodynamics_2d.hpp
  \brief Various manipulations of hydrodynamic variables
  \author Almog Yalinewich
*/

#ifndef HYDRODYNAMICS_2D_HPP
#define HYDRODYNAMICS_2D_HPP 1

#include "spatial_distribution2d.hpp"
#include "../common/equation_of_state.hpp"
#include "point_motion.hpp"
#include "spatial_reconstruction.hpp"
#include "SourceTerm.hpp"
#include "CustomEvolution.hpp"
#include "../../misc/utils.hpp"
#include "../common/hydrodynamics.hpp"
#include "../../misc/universal_error.hpp"
#include "../../tessellation/ConvexHull.hpp"
#include "../../tessellation/polygon_overlap_area.hpp"
#include <boost/scoped_ptr.hpp>
#include "../../tessellation/EdgeLengthCorrect.hpp"
#include "../../mpi/mpi_macro.hpp"
#include "../../mpi/ProcessorUpdate.hpp"
#include "physical_geometry.hpp"

/*! \brief Initialize computational cells
  \param density Density distribution
  \param pressure Pressure distribution
  \param xvelocity Distribution of the x component of the velocity
  \param yvelocity Distribution of the y component of the velocity
  \param eos Equation of state
  \param tess Tessellation
  \param CMvalue Determines whether to evaluate the spatial distributions in the mesh generating points or the center of mass
  \return List of primitive variables
*/
vector<Primitive> InitialiseCells
(SpatialDistribution const& density,
 SpatialDistribution const& pressure,
 SpatialDistribution const& xvelocity,
 SpatialDistribution const& yvelocity,
 EquationOfState const& eos,
 Tessellation const& tess,bool CMvalue=true);

/*! \brief Calculates the intensive conserved variables
  \param cells Hydrodynamical cells
  \return List of conserved variables
*/
vector<Conserved> CalcConservedIntensive
(vector<Primitive> const& cells);

/*! \brief Calculates the extensive conserved variables
  \param cons_int Conserved intensive variables
  \param tess Tessellation
  \param pg Physical geometry
  \return List of conserved variables
*/
vector<Conserved> CalcConservedExtensive
(const vector<Conserved>& cons_int,
 const Tessellation& tess,
 const PhysicalGeometry& pg);

/*! \brief Calculates the velocities of the mesh generating points
  \param tessellation Tessellation
  \param cells Computational cells
  \param pointmotion Point motion function
  \param pointvelocity List of velocities
  \param time Time
	\param cevolve The custom evolution of the cells
*/
void CalcPointVelocities(Tessellation const& tessellation,
			 vector<Primitive> const& cells,
			 PointMotion& pointmotion,
			 vector<Vector2D>& pointvelocity,double time,
			vector<CustomEvolution*> & cevolve);

/*! \brief Calculates the time step for a cell
  \param cell Computational cell
  \param width Cell width
  \param face_velocites Velocities of the edges of the cell
  \return Time step
*/
double TimeStepForCell(Primitive const& cell,
		       double width,vector<Vector2D> const& face_velocites);

/*! \brief Calculates the time step
  \param tessellation Tessellation
  \param cells Hydrodynamic fluid elements
  \param facevelocity Velocities of the interfaces between cells
  \param hbc Hydrodynamic boundary condition
  \param time Time
  \param evolve Custom evolution
  \return Time step
*/
double CalcTimeStep
(Tessellation const& tessellation,
 vector<Primitive> const& cells,
 vector<Vector2D> const& facevelocity,
 HydroBoundaryConditions const& hbc,
 double time,
 vector<CustomEvolution*> const& evolve=vector<CustomEvolution*>());

/*! \brief Updates the extensive conserved variables
  \param tessellation Tessellation
  \param fluxes Hydrodynamic fluxes
  \param dt Time step
  \param conserved_extensive Extensive conserved variables
  \param boundaryconditions Hydrodynamic boundary condition
  \param lengthes The corrected lengthes of the edges
*/
void UpdateConservedExtensive
(Tessellation const& tessellation,
 vector<Conserved> const& fluxes,
 double dt,
 vector<Conserved>& conserved_extensive,
 HydroBoundaryConditions const& boundaryconditions,
 vector<double> const& lengthes);

/*! \brief Move mesh points
  \param pointvelocity Velocities of all mesh points
  \param dt Time step
  \param tessellation Tessellation
  \param oldpoints Possible input for the location of the old points
*/
void MoveMeshPoints(vector<Vector2D> const& pointvelocity,
		    double dt, Tessellation& tessellation,
			vector<Vector2D> oldpoints=vector<Vector2D> ());
/*! \brief Move mesh points
  \param pointvelocity Velocities of all mesh points
  \param dt Time step
  \param tessellation The local tessellation
  \param vproc The processor tessellation
  \param oldpoints Possible input for the location of the old points
*/
void MoveMeshPoints(vector<Vector2D> const& pointvelocity,
	double dt, Tessellation& tessellation,Tessellation const& vproc,
	vector<Vector2D> oldpoints=vector<Vector2D> ());

/*! \brief Calculates the intensive conserved variables
  \param tess Tessellation
  \param extensive Extensive conserved variables
  \param pg Physical geometry
  \return List of intensive conserved variables
 */
vector<Conserved> calc_conserved_intensive
(const Tessellation& tess,
 const vector<Conserved>& extensive,
 const PhysicalGeometry& pg);

/*! \brief Updates the intensive conserved variables
  \param tessellation Tessellation
  \param conservedextensive Extensive conserved variables
  \param conservedintensive Intensive conserved variables
*/
void UpdateConservedIntensive(Tessellation const& tessellation,
			      vector<Conserved> const& conservedextensive,
			      vector<Conserved>& conservedintensive);

/*! \brief Updates the primitive variables
  \param conservedintensive Intensive conserved variables
  \param eos Equation of state
  \param cells Fluid elements
  \param CellsEvolve Custom evolution
  \param old_cells Fluid element in the previous half step
  \param densityfloor Toggle density floor protection
  \param densitymin Minimum density
  \param pressuremin Minimum pressure
  \param tess Tessellation
  \param time Time
  \param extensivetracers Extensive tracers
  \todo Pass old_cells as const and encapsulate densitymin parameters
*/
void UpdatePrimitives(vector<Conserved> const& conservedintensive,
		      EquationOfState const& eos,vector<Primitive>& cells,
		      vector<CustomEvolution*> const& CellsEvolve,
		      vector<Primitive> &old_cells,bool densityfloor,
		      double densitymin,double pressuremin,
		      Tessellation const& tess,double time,
		      vector<vector<double> > const& extensivetracers);

/*! \brief Calculates the fluxes
  \param tessellation Tessellation
  \param cells List of primitive variables
  \param dt Time step
  \param time Simulation time
  \param interpolation Interpolation scheme
  \param facevelocity Velocity of edges
  \param boundaryconditions Hydrodynamic boundary conditions
  \param rs Riemann solver
  \param CellsEvolve Custom evolution schemes
  \param cem Custom evolution manager
  \param tracers Tracers
  \return List of fluxes at each edge
 */
vector<Conserved> calc_fluxes
(Tessellation const& tessellation,
 vector<Primitive> const& cells,
 double dt,
 double time,
 SpatialReconstruction& interpolation,
 vector<Vector2D> const& facevelocity,
 HydroBoundaryConditions const& boundaryconditions,
 RiemannSolver const& rs,
 vector<CustomEvolution*> const& CellsEvolve,
 CustomEvolutionManager const& cem,
 vector<vector<double> > const& tracers);

/*! \brief Calculates the flux in the bulk of the fluid
  \param normaldir A unit vector normal to the interface
  \param paraldir A unit vector parallel to the interface
  \param left Primitive variables on the left side of the interface
  \param right Primitive variables on the right side of the interface
  \param edge_velocity Velocity of the interface
  \param rs Riemann solver
  \return Flux
 */
Conserved FluxInBulk(Vector2D const& normaldir,
		     Vector2D const& paraldir,
		     Primitive const& left,
		     Primitive const& right,
		     Vector2D const& edge_velocity,
		     RiemannSolver const& rs);

/*! \brief Adds force contribution to the extensive conserved variables
  \param tess Tessellation
  \param pg Physical geometry
  \param cells Fluid elements
  \param force External force
  \param t Time
  \param dt Time step
  \param conserved_extensive Extensive conserved variables
  \param hbc Hydrodynamic boundary conditions
  \param fluxes Hydrodynamic fluxes
  \param point_velocity Velocities of the mesh generating points
  \param g Output
  \param coldflows_flag Toggles cold flows
  \param tracers Tracers
  \param lengthes The corrected length of the edges
 */
void ExternalForceContribution(Tessellation const& tess,
			       const PhysicalGeometry& pg,
			       vector<Primitive> const& cells,
			       SourceTerm& force,
			       double t,
			       double dt,
			       vector<Conserved>& conserved_extensive,
			       HydroBoundaryConditions const& hbc,
			       vector<Conserved> const& fluxes,
			       vector<Vector2D> const& point_velocity,
			       vector<double> &g,
			       bool coldflows_flag,vector<vector<double> > &tracers,
				   vector<double> const& lengthes);

/*! \brief Second order time advance
  \param tess Tessellation
  \param proctess The tessellation of the processors (if no mpi pass the regular tesselation)
  \param cells Fluid elements,
  \param point_motion Point motion scheme
  \param hbc Hydrodynamic boundary conditions
  \param interpolation Spatial reconstruction
  \param rs Riemann solver
  \param eos Equation of state
  \param force External source term
  \param time Time
  \param cfl Courant Friedrich Lewy number
  \param endtime Final time for the simulation
  \param tracers Tracers
  \param dt_external Extrnal time step
  \param custom_evolution_indices The indices of the customevolution
  \param custom_evolution_manager Class that translates indices to class pointers
  \param pg Physical geometry
  \param procupdate Scheme for updating the positions of the processes
  \param traceflag Determines whether tracers should be updated
  \param coldflows_flag Determines whether cold flows should be used
  \param as Described in the Arepo paper
  \param bs Described in the Arepo paper
  \param densityfloor Determines whether densityfloor should be used
  \param densitymin Minimum density
  \param pressuremin Minimum pressure
  \param EntropyCalc Determines whether entropy should be calculated
  \return Time step
  \todo document and encapsulate parameters as, bs
 */
double TimeAdvance2mid
(Tessellation& tess,
#ifdef RICH_MPI
 Tessellation& proctess,
#endif
 vector<Primitive> &cells,
 PointMotion& point_motion,
 HydroBoundaryConditions const& hbc,
 SpatialReconstruction& interpolation,
 RiemannSolver const& rs,
 EquationOfState const& eos,
 SourceTerm& force,
 double time,
 double cfl,
 double endtime,
 vector<vector<double> >& tracers,
 double dt_external,
 vector<size_t>& custom_evolution_indices,
 const CustomEvolutionManager& custom_evolution_manager,
 const PhysicalGeometry& pg,
 #ifdef RICH_MPI
 ProcessorUpdate *procupdate,
 #endif
 bool traceflag=false,
 bool coldflows_flag=false,double as=0.01,
 double bs=0.01,
 bool densityfloor=false,
 double densitymin=0.01,
 double pressuremin=0.01,
 bool EntropyCalc=false);

/*!
  \brief Essential data needed to advance a numerical simulation to the next time step
*/
/*! \brief Calculates the velocities of the mesh generating points
  \param tess Tessellation
  \param cells Fluid elements
  \param point_motion Point motion scheme
  \param time Time
\param cevolve The custom evolution of the cells
  \return List of the velocities of the mesh generating points
 */
vector<Vector2D> calc_point_velocities
(Tessellation const& tess,
 vector<Primitive>const& cells,
 PointMotion& point_motion,
 double time,vector<CustomEvolution*> & cevolve);

/*! \brief Returns the position of all mesh generating points
  \param tess Tessellation
  \return Position of all mesh generating points
 */
vector<Vector2D> get_all_mesh_points
(Tessellation const& tess);

/*! \brief Changes the energy and sound speed so they would satisfy the equation of state
  \param vp Primitive variables
  \param eos Equation of state
  \return Primitive variable with corrected sound speed and energy
 */
vector<Primitive> make_eos_consistent
(vector<Primitive> const& vp,
 EquationOfState const& eos);

/*! \brief Calculates the change in the tracer
  \param old_trace Tracers in the beginning of the time step
  \param cells Fluid elements
  \param tess Tessellation
  \param fluxes Hydrodynamic fluxes
  \param dt Time step
  \param bc Hydrodynamic boundary conditions
  \param interp Spatial reconstruction
  \param time Time
  \param cellsevolve Custom evolution
  \param cem Custom evolution manager
  \param edge_velocities Velocities of the edges between cells
  \param lengthes The corrected lengthes of the edges
  \return List of tracers
 */
vector<vector<double> > CalcTraceChange
(vector<vector<double> > const& old_trace,
 vector<Primitive> const& cells,
 Tessellation const& tess,vector<Conserved> const& fluxes,double dt,
 HydroBoundaryConditions const& bc,
 SpatialReconstruction const& interp,
 double time,vector<CustomEvolution*> const& cellsevolve,
 CustomEvolutionManager const& cem,
 vector<Vector2D> const& edge_velocities,vector<double> const& lengthes);

/*! \brief Calculates the maximum kinetic energy
  \param tess Tessellation
  \param cells Fluid elements
  \param customevolve Custom evolution
  \return List of kinetic energies
 */
vector<double> GetMaxKineticEnergy
(Tessellation const& tess,vector<Primitive> const& cells,
 vector<CustomEvolution*> const& customevolve);

/*! \brief Returns the energies due to external potentials
  \param tess Tessellation
  \param g TBA
  \return TBA
  \todo Add documentation
 */
vector<double> GetForceEnergy(Tessellation const& tess,
			      vector<double> const& g);

/*! \brief Cold flows pressure fix
  \param intensive Intensive conserved variables
  \param entropy Entropy
  \param eos Equation of state
  \param Ek See Arepo paper
  \param Ef See Arepo paper
  \param as See Arepo paper
  \param bs See Arepo paper
  \param customevolve Custom evolution
  \param tess Tessellation
  \param shockedcells Flag whether the cell is shocked or not
  \param densityfloor Is the densityfllor on or off
 */
void FixPressure
(vector<Conserved> &intensive,vector<vector<double> > const& entropy,
 EquationOfState const& eos,vector<double> const& Ek,
 vector<double> const& Ef,double as,double bs,vector<CustomEvolution*>
 const& customevolve,Tessellation const& tess,/*vector<Conserved>
 &extensive,*/vector<char> const& shockedcells,bool densityfloor);

/*! \brief Returns true is a cell is near a boundary
  \param index Cell index
  \param tess Tessellation
  \param customevolve Custom evolusion
  \return True if a cell is near a boundary
  \todo Remove custom evolution from the signature
 */
bool NearBoundary(int index,Tessellation const& tess,
		  vector<CustomEvolution*> const& customevolve);

/*! \brief Calculates extensive tracers
  \param intensive_tracer Intensive tracers
  \param tess Tessellation
  \param cells List of primitive variables
  \param pg Physical geometry
  \return List of extensive tracers
 */
vector<vector<double> > calc_extensive_tracer
(const vector<vector<double> > & intensive_tracer,
 const Tessellation& tess,
 const vector<Primitive>& cells,
 const PhysicalGeometry& pg);

/*! \brief Calculates the extensive tracer
  \param tracer Intensive tracer
  \param tess Tessellation
  \param cells Fluid elements
  \param result List of extensive tracers
 */
void MakeTracerExtensive
(vector<vector<double> > const&tracer,
 Tessellation const& tess,vector<Primitive> const& cells,
 vector<vector<double> > &result);

/*! \brief Calculates the intensive tracers
  \param tracer Output
  \param tracer_extensive Extensive tracers
  \param tess Tessellation
  \param cells Fluid elements
  \param pg Physical geometry
 */
void MakeTracerIntensive
(vector<vector<double> > &tracer,vector<vector<double> >
 const& tracer_extensive,Tessellation const& tess,
 vector<Primitive> const& cells, const PhysicalGeometry& pg);

/*! \brief A more efficient version of update_extensive_tracers
  \param extensive_tracers Extensive tracers
  \param tracers Intensive tracers
  \param cells List of primitive variables
  \param tess Tessellation
  \param fluxes Hydrodynamic fluxes
  \param time Simulation time
  \param dt Time step
  \param hbc Hydrodynamic variabels
  \param interp Interpolation scheme
  \param ce Custom evolutions
  \param cem Custom evolution manager
  \param fv Face velocities
  \param lengths Edge lengths
 */
void really_update_extensive_tracers
(vector<vector<double> >& extensive_tracers,
 const vector<vector<double> >& tracers,
 const vector<Primitive>& cells,
 const Tessellation& tess,
 const vector<Conserved>& fluxes,
 double time, double dt,
 const HydroBoundaryConditions& hbc,
 const SpatialReconstruction& interp,
 const vector<CustomEvolution*>& ce,
 const CustomEvolutionManager& cem,
 const vector<Vector2D>& fv,
 const vector<double>& lengths);

/*! \brief Updates the extensive tracers
  \param tracerextensive Extensive tracers
  \param tracerchange Change in extensive tracers
  \param CellsEvolve Custom evolution
  \param cells Fluid elements
  \param tess Tessellation
  \param time Time
 */
void UpdateTracerExtensive
(vector<vector<double> > &tracerextensive,
 vector<vector<double> > const& tracerchange,vector<CustomEvolution*> const&
 CellsEvolve,vector<Primitive> const& cells,Tessellation const& tess,double time);

/*! \brief Resets the primitive variables based on a tracer threshold
  \param alpha The tracer threshold
  \param originalD The original density distribution.
  \param originalP The original pressure distribution.
  \param originalVx The original x velocity distribution.
  \param originalVy The original y velocity distribution.
  \param cells The primitive cells
  \param tess The tessellation
  \param tracer The tracer field
  \param tracerindex The index in the tracer to consider
  \param eos The equation of state
  \param cevolve The custom evolution of cells
  \param coldflows Toggles cold flows correction
*/
void TracerResetCalc(double alpha,SpatialDistribution const& originalD,
		     SpatialDistribution const& originalP,SpatialDistribution const& originalVx,
			 SpatialDistribution const& originalVy, vector<SpatialDistribution const*> const& ,
			 vector<Primitive> &cells,Tessellation const& tess,vector<vector<double> > &tracer,
		     int tracerindex,EquationOfState const& eos,vector<CustomEvolution*>
			const& cevolve,bool coldflows);

/*! \brief Makes a list of points to remove
  \param tess Tessellation
  \param point TBA
  \param R Radius?
  \param PointToRemove output
  \param Inner TBA
  \todo Add documentation
 */
void GetPointToRemove(Tessellation const& tess,Vector2D const& point,
		      double R,vector<int> & PointToRemove,int Inner);

/*! \brief Returns true if a computational cell is shocked
  \param tess Tessellation
  \param index Cell index
  \param cells Fluid elements
  \param hbc Hydrodynamic boundary conditions
  \param time Time
  \return True if a cell is shocked
 */
bool IsShockedCell
(Tessellation const& tess,int index,
 vector<Primitive> const& cells,HydroBoundaryConditions const& hbc,
 double time);

/*! \brief Calculates time step for a cell boundary
  \param cell Fluid element
  \param cells Fluid elements
  \param width Cell width
  \param face_velocities Velocities of the cell interfaces
  \param tess Tessellation
  \param hbc Hydrodynamic boundary conditions
  \param index Cells index
  \param time Time
  \return Time step
  \todo Check for redundancy (can't cell be reproduced from cells and index?)
 */
double TimeStepForCellBoundary
(Primitive const& cell,
 vector<Primitive> const& cells,double width,
 vector<Vector2D> const& face_velocities,Tessellation const& tess,
 HydroBoundaryConditions const& hbc,int index,double time);

/*!
\brief Convert indeces of customevolution to pointer to classes
\param cem The class that relates indeces to customevolution
\param indices The indices of the customevolution
\return The pointers to the customevolution classes
*/
vector<CustomEvolution*> convert_indices_to_custom_evolution(const CustomEvolutionManager& cem,
	const vector<size_t>& indices);

/*! \brief Applies a correction to the extensive variables due to the change in volume during time step.
  \details This method calculates the change in extensive by calculating the volume swept by an edge and multiplying it by the intensive variables of the respective cell.
  \param extensive Extensive variables
  \param intensive Intensive variables
  \param tessold Old tessellation
  \param tessnew New tessellation
  \param facevelocity Face velocity
  \param dt Time step
  \param pointvelocity Velocities of mesh generating points
 */
void FixAdvection(vector<Conserved>& extensive,
		  vector<Conserved> const& intensive,
		  Tessellation const& tessold,
		  Tessellation const& tessnew,
		  vector<Vector2D> const& facevelocity,
		  double dt,vector<Vector2D> const& pointvelocity);

/*! \brief Determines the time step
  \param hydro_time_step Time step derived from hydrodynamics
  \param external_dt Time step suggested by user
  \param current_time Current simulation time
  \param end_time Termination time
  \return Time step
 */
double determine_time_step(double hydro_time_step,
			   double external_dt,
			   double current_time,
			   double end_time);

#endif // HYDRODYNAMICS_2D_HPP
