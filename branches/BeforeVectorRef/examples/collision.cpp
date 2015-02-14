#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/misc/int2str.hpp"
#include "collide_density.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"

class ProgressReport: public DiagnosticFunction
{
public:

  ProgressReport(int skip):
    skip_(skip) {}

  void operator()(const hdsim& sim)
  {
    if(sim.GetCycle()%skip_==0)
      cout << "Sim time is " << sim.GetTime() << " Step number "
	   << sim.GetCycle() << endl;
  }

private:
  const int skip_;
};

int main(void)
{
  // Set up the initial grid points
  double blob_radius=0.2;
  int npoints_blob=1000;
  double blob1_xc=-0.5;
  double blob2_xc=0.5;
  double blob_yc=0;

  // Set up the boundary type for the points
  SquareBox outer(Vector2D(-1,-1), Vector2D(1,1));

  vector<Vector2D> background = 
    cartesian_mesh(50,50, 
		   outer.getBoundary().first, 
		   outer.getBoundary().second);
  vector<Vector2D> blob1=CirclePointsRmax(npoints_blob,blob_radius*0.1,
					  blob_radius,blob1_xc,blob_yc);
  vector<Vector2D> blob2=CirclePointsRmax(npoints_blob,blob_radius*0.1,
					  blob_radius,blob2_xc,blob_yc);
  vector<Vector2D> InitPoints=background;
  InitPoints.insert(InitPoints.end(),blob1.begin(),blob1.end());
  InitPoints.insert(InitPoints.end(),blob2.begin(),blob2.end());

  // Set up the tessellation
  VoronoiMesh tess;

  // Set up the Riemann solver
  Hllc rs;

  // Set the hydro boundary conditions
  RigidWallHydro hbc(rs);

  // Set up the equation of state
  IdealGas eos(5./3.);

  // Set up the point motion scheme
  Lagrangian l_motion;
  RoundCells pointmotion(l_motion,hbc);

  // Set up the interpolation
  LinearGaussConsistent interpolation(eos,outer,hbc);

  // Set up the initial Hydro
  CollideDensity density(10,10,0.1,blob1_xc,blob_yc,
			 blob2_xc,blob_yc,blob_radius,blob_radius);
  Uniform2D pressure(1);
  CollideDensity xvelocity(1,-1,0,blob1_xc,
			   blob_yc,blob2_xc,blob_yc,blob_radius,blob_radius);
  Uniform2D yvelocity(0);

  // Set up the external source term
  ZeroForce force;

  // Set up the simulation
  hdsim sim(InitPoints,tess,interpolation,density,pressure,xvelocity,
	    yvelocity,eos,rs,pointmotion,force,outer,hbc);

  // Choose the Courant number
  sim.SetCfl(0.7);

  // How long shall we run the simulation?
  double tend=1;
  sim.SetEndTime(tend);

  // Custom output criteria
  SafeTimeTermination term_cond(tend,1e6);
  ProgressReport diag1(25);
  ConsecutiveSnapshots diag2(0.1);
  MultipleDiagnostics diag;
  diag.diag_list.push_back(&diag1);
  diag.diag_list.push_back(&diag2);
  simulation2d::main_loop(sim,
			  term_cond,
			  &hdsim::TimeAdvance2Mid,
			  &diag);

  // Done running the simulation, output the data
  write_snapshot_to_hdf5(sim,"final.h5");

  // We are done!!
  cout<<"Finished running the simulation"<<endl;

  return 0;
}
