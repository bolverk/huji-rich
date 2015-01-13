#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylinderical_geometry.hpp"
#include "source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/tessellation/shape_2d.hpp"
#include "source/newtonian/test_2d/piecewise.hpp"
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
  // Set up the boundary type for the points
  SquareBox outer(Vector2D(0,0),
		  Vector2D(1,1));


  // Set up the initial grid points
  vector<Vector2D> InitPoints = 
    cartesian_mesh(50,50,
		   outer.getBoundary().first,
		   outer.getBoundary().second);

  // Set up the tessellation
  VoronoiMesh tess;

  // Set up the Riemann solver
  Hllc rs;

  // Set the hydro boundary conditions
  RigidWallHydro hbc(rs);

  // Set up the equation of state
  IdealGas eos(5./3.);

  // Set up the point motion scheme
  Lagrangian raw_pointmotion;
  RoundCells pointmotion(raw_pointmotion, hbc);  

  // Set up the interpolation
  LinearGaussConsistent interpolation(eos,outer,hbc);

  // Set up the initial Hydro
  Uniform2D density(1);
  Circle hot_spot(Vector2D(0,0),0.3);
  Uniform2D low_pressure(1);
  Uniform2D high_pressure(100);
  Piecewise pressure(hot_spot,
		     high_pressure,
		     low_pressure);
  Uniform2D xvelocity(0);
  Uniform2D yvelocity(0);

  // Set up the external source term
  CylindericalGeometry force(Vector2D(0,0),Vector2D(0,1));

  // Set up the simulation
  hdsim sim(InitPoints,tess,interpolation,density,pressure,xvelocity,
	    yvelocity,eos,rs,pointmotion,force,outer,hbc);

  // Add the tracer for following the mass
  Uniform2D tracer_inside(1);
  Uniform2D tracer_outside(0);
  Piecewise tracer(hot_spot,
		   tracer_inside,
		   tracer_outside);
  sim.addTracer(tracer);

  // Choose the Courant number
  sim.SetCfl(0.7);

  // How long shall we run the simulation?
  const double tend=0.05;
  sim.SetEndTime(tend);

  // Main loop
  SafeTimeTermination term_cond(tend,1e6);
  ProgressReport diag1(25);
  ConsecutiveSnapshots diag2(0.01);
  MultipleDiagnostics diag;
  diag.diag_list.push_back(&diag1);
  diag.diag_list.push_back(&diag2);
  simulation2d::main_loop(sim,
			  term_cond,
			  &hdsim::TimeAdvance2Mid,
			  &diag);

  // Done running the simulation, output the data
  write_snapshot_to_hdf5(sim,"final.h5");
	
  // Announce 
  cout<<"Finished running the simulation"<<endl;

  return 0;
}
