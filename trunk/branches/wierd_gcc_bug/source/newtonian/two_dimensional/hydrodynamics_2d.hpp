#ifndef HYDRODYNAMICS_2D_HPP
#define HYDRODYNAMICS_2D_HPP 1

#include <vector>
#include "../common/hydrodynamic_variables.hpp"
#include "spatial_distribution2d.hpp"
#include "../common/equation_of_state.hpp"
#include "point_motion.hpp"
#include "spatial_reconstruction.hpp"
#include "SourceTerm.hpp"
#include "CustomEvolution.hpp"
#include "scalar_interpolation.hpp"
#include "../../misc/utils.hpp"

using namespace std;

vector<Primitive> InitialiseCells
	(SpatialDistribution const& density,
	SpatialDistribution const& pressure,
	SpatialDistribution const& xvelocity,
	SpatialDistribution const& yvelocity,
	EquationOfState const& eos,
	Tessellation const* tess,bool CMvalue=true);

vector<Conserved> CalcConservedIntensive
	(vector<Primitive> const& cells);

vector<Conserved> CalcConservedExtensive
	(vector<Conserved> const& cons_int,
	Tessellation const* tess);

void CalcPointVelocities(Tessellation const* tessellation,
	vector<Primitive> const& cells,
	PointMotion *pointmotion,
	vector<Vector2D>& pointvelocity,double time);

double TimeStepForCell(Primitive const& cell,
	double width,vector<Vector2D> const& face_velocites);

double CalcTimeStep(Tessellation const* tessellation,
	vector<Primitive> const& cells,
	vector<Vector2D> const& facevelocity,
	HydroBoundaryConditions const* hbc,double time);

void UpdateConservedExtensive
	(Tessellation const* tessellation,
	vector<Conserved> const& fluxes,
	double dt, 
	vector<Conserved>& conserved_extensive,
	HydroBoundaryConditions const* boundaryconditions);

void MoveMeshPoints(vector<Vector2D> const& pointvelocity,
	double dt, Tessellation* tessellation);

void UpdateConservedIntensive(Tessellation const* tessellation,
	vector<Conserved> const& conservedextensive,
	vector<Conserved>& conservedintensive);

/*void UpdateConservedIntensive(Tessellation const* tessellation,
vector<Conserved> const& conservedextensive,
vector<Conserved>& conservedintensive);*/

void UpdatePrimitives(vector<Conserved> const& conservedintensive,
	EquationOfState const& eos,vector<Primitive>& cells,
	vector<CustomEvolution*> const& CellsEvolve,
	vector<Primitive> &old_cells,bool densityfloor=false,
	double densitymin=0.01,double pressuremin=0.01);

void CalcFluxes
	(Tessellation const* tessellation,
	vector<Primitive> const& cells,
	double dt,
	double time,
	SpatialReconstruction* interpolation,
	vector<Vector2D> const& facevelocity,
	HydroBoundaryConditions const* boundaryconditions,
	RiemannSolver const& rs,
	vector<Conserved>& fluxes,
	vector<CustomEvolution*> const& CellsEvolve,
	vector<vector<double> > const& tracers);

Conserved FluxInBulk(Vector2D const& normaldir,
	Vector2D const& paraldir,
	Primitive const& left,
	Primitive const& right,
	Vector2D const& edge_velocity,
	RiemannSolver const& rs);

void ExternalForceContribution(Tessellation const* tess,
	vector<Primitive> const& cells,
	SourceTerm *force,
	double t,
	double dt,
	vector<Conserved>& conserved_extensive,
	HydroBoundaryConditions const* hbc,
	vector<Conserved> const& fluxes,
	vector<Vector2D> const& point_velocity,
	vector<double> &g, 
	bool coldflows_flag,vector<vector<double> > &tracers);

double TimeAdvance2mid(Tessellation* tess,
	vector<Primitive> &cells,PointMotion *point_motion,
	HydroBoundaryConditions const* hbc,
	SpatialReconstruction* interpolation,
	RiemannSolver const& rs,EquationOfState const& eos,
	SourceTerm* force,double time,double cfl,double endtime,
	vector<CustomEvolution*> const& CellsEvolve,vector<vector<double> >
	&tracers,double dt_external,
	bool traceflag=false,bool coldflows_flag=false,double as=0.01,
	double bs=0.01,bool densityfloor=false,double densitymin=0.01,
	double pressuremin=0.01,bool EntropyCalc=false);

/*! 
\brief Essential data needed to advance a numerical simulation to the next time step
*/
class HydroSnapshot
{
public:
	/*! \brief Class constructor
	\param mesh_points The mesh points
	\param cells The primitive cells
	*/
	HydroSnapshot(vector<Vector2D> const& mesh_points,
		vector<Primitive> const& cells);
	//! \brief Default constructor
	HydroSnapshot();
	//! \brief The mesh points
	vector<Vector2D> mesh_points;
	//! \brief The primitive cells
	vector<Primitive> cells;  
};

vector<Vector2D> calc_point_velocities
	(Tessellation const* tess,
	vector<Primitive>const& cells,
	PointMotion *point_motion,
	double time);

vector<Vector2D> get_all_mesh_points
	(Tessellation const* tess);

vector<Primitive> make_eos_consistent
	(vector<Primitive> const& vp,
	EquationOfState const& eos);

vector<vector<double> > CalcTraceChange
(vector<vector<double> > const& old_trace,vector<Primitive> const& cells,
 Tessellation const* tess,vector<Conserved> const& fluxes,double dt,
 HydroBoundaryConditions const* bc,
 SpatialReconstruction const* interp,
 double time);

vector<double> GetMaxKineticEnergy(Tessellation const* tess,vector<Primitive> const&
	cells,vector<CustomEvolution*> const& customevolve);

vector<double> GetForceEnergy(Tessellation const* tess,
	vector<double> const& g);

void FixPressure(vector<Conserved> &intensive,vector<vector<double> > const& entropy,
	EquationOfState const& eos,vector<double> const& Ek,
	vector<double> const& Ef,double as,double bs,vector<CustomEvolution*>
	const& customevolve,Tessellation const* tess,vector<Conserved>
	&extensive,vector<Primitive> const& cells,
	HydroBoundaryConditions const* hbc,double time);

bool NearBoundary(int index,Tessellation const* tess,
	vector<CustomEvolution*> const& customevolve);

void MakeTracerExtensive(vector<vector<double> > const&tracer,
	Tessellation const* tess,vector<Primitive> const& cells,vector<vector<double> >
	&result);

void MakeTracerIntensive(vector<vector<double> > &tracer,vector<vector<double> >
	const& tracer_extensive,Tessellation const* tess,vector<Primitive> const& cells);

void UpdateTracerExtensive(vector<vector<double> > &tracerextensive,
	vector<vector<double> > const& tracerchange);
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
\param InnerNum The number of cells to always reset
*/
void TracerResetCalc(double alpha,SpatialDistribution const& originalD,
	SpatialDistribution const& originalP,SpatialDistribution const& originalVx,
	SpatialDistribution const& originalVy, vector<Primitive> &cells,
	Tessellation const* tess,vector<vector<double> > &tracer,
	int tracerindex,EquationOfState const& eos,int InnerNum);

void GetPointToRemove(Tessellation const* tess,Vector2D const& point,
	double R,vector<int> & PointToRemove,int Inner);

bool IsShockedCell(Tessellation const* tess,int index,
	vector<Primitive> const& cells,HydroBoundaryConditions const* hbc,
	double time);

double TimeStepForCellBoundary(Primitive const& cell,
	vector<Primitive> const& cells,double width,
	vector<Vector2D> const& face_velocites,Tessellation const* tess,
	HydroBoundaryConditions const* hbc,int index,double time);

#endif // HYDRODYNAMICS_2D_HPP
