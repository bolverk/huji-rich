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
#include "../../misc/utils.hpp"
#include "../../misc/lazy_list.hpp"
#include "../common/hydrodynamics.hpp"
#include "../../misc/universal_error.hpp"
#include "../../tessellation/ConvexHull.hpp"
#include "../../tessellation/polygon_overlap_area.hpp"
#include <boost/scoped_ptr.hpp>
#include "../../tessellation/EdgeLengthCorrect.hpp"
#include "physical_geometry.hpp"
#include "time_step_function.hpp"
#include "cache_data.hpp"
#include "../common/riemann_solver.hpp"

/*! \brief Rotates primitive variables to align with edge
  \param n Normal directions
  \param p Parallel direction
  \param cell Primitive variables
  \return Rotated cell
 */
Primitive RotatePrimitive(const Vector2D& n,
	const Vector2D& p,
	const Primitive& cell);

/*! \brief Rotates flux from the edge frame back to the lab frame
  \param c Flux
  \param n Normal direction
  \param p Parallel direction
  \return Rotated flux
 */
Conserved RotateFluxBack(const Conserved& c,
	const Vector2D& n,
	const Vector2D& p);

/*! \brief Given an edge and an index of one neighbor, returns the index of another neighbor
  \param edge Voronoi edge
  \param index Index of one neighbor
  \return Index of the other neighbor
 */
int get_other_index(const Edge& edge,
	const int index);

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
	Tessellation const& tess, bool CMvalue = true);

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

/*! \brief Calculates the time step for a cell
  \param cell Computational cell
  \param width Cell width
  \param face_velocites Velocities of the edges of the cell
  \return Time step
*/
double TimeStepForCell(Primitive const& cell,
	double width, vector<Vector2D> const& face_velocites);

/*! \brief Move mesh points
  \param pointvelocity Velocities of all mesh points
  \param dt Time step
  \param tessellation Tessellation
  \param oldpoints Possible input for the location of the old points
  \param reorder Should we reorder the points
  \return The indeces of the hilbert order of the points
*/
vector<int> MoveMeshPoints(vector<Vector2D> const& pointvelocity,
	double dt, Tessellation& tessellation, bool reorder,
	vector<Vector2D> oldpoints = vector<Vector2D>());
/*! \brief Move mesh points
  \param pointvelocity Velocities of all mesh points
  \param dt Time step
  \param tessellation The local tessellation
  \param vproc The processor tessellation
  \param oldpoints Possible input for the location of the old points
  \param reorder Should we reorder the points
  \return The indeces of the hilbert order of the points
*/
vector<int> MoveMeshPoints(vector<Vector2D> const& pointvelocity,
	double dt, Tessellation& tessellation, Tessellation const& vproc, bool reorder,
	vector<Vector2D> oldpoints = vector<Vector2D>());

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
  \param cd Cache data
  \param cells Computational cells
  \param fluxes Fluxes
  \param point_velocities Velocities of the mesh generating points
  \param source Source term
  \param t Time
  \param dt Time step
  \param extensives Extensive variables
  \param tracerstickernames The names of the tracers and stickers
 */
void ExternalForceContribution
(const Tessellation& tess,
	const PhysicalGeometry& pg,
	const CacheData& cd,
	const vector<ComputationalCell>& cells,
	const vector<Extensive>& fluxes,
	const vector<Vector2D>& point_velocities,
	const SourceTerm& source,
	double t,
	double dt,
	vector<Extensive>& extensives,
	TracerStickerNames const& tracerstickernames);

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

/*! \brief Returns the energies due to external potentials
  \param tess Tessellation
  \param g TBA
  \return TBA
  \todo Add documentation
 */
vector<double> GetForceEnergy(Tessellation const& tess,
	vector<double> const& g);

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
	Tessellation const& tess, vector<Primitive> const& cells,
	vector<vector<double> > &result);

/*! \brief Makes a list of points to remove
  \param tess Tessellation
  \param point TBA
  \param R Radius?
  \param PointToRemove output
  \param Inner TBA
  \todo Add documentation
 */
void GetPointToRemove(Tessellation const& tess, Vector2D const& point,
	double R, vector<int> & PointToRemove, int Inner);

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
	double dt, vector<Vector2D> const& pointvelocity);

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
