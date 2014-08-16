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
#include "source/tessellation/shape_2d.hpp"
#include "source/newtonian/test_2d/piecewise.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"

namespace {
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

  void HeatCells(hdsim &sim,double minDensity,
		 double newDensity,double newPressure, 
		 EquationOfState const& eos)
  {
    vector<Primitive>& cells=sim.GetAllCells();
    int n=(int)cells.size();
    for(int i=0;i<n;++i)
      {
	if(cells[i].Density<minDensity)
	  {
	    cells[i].Density=newDensity;
	    cells[i].Pressure=newPressure;
	    cells[i].SoundSpeed=eos.dp2c(newDensity,newPressure);
	    cells[i].Energy=eos.dp2e(newDensity,newPressure);
	  }
      }
  }

  class CellHeater: public Manipulate
  {
  public:

    CellHeater(double min_density,
	       double new_density,
	       double new_pressure,
	       const EquationOfState& eos):
      min_density_(min_density),
      new_density_(new_density),
      new_pressure_(new_pressure),
      eos_(eos) {}

    void operator()(hdsim& sim)
    {
      HeatCells(sim,min_density_,new_density_,new_pressure_,eos_);
    }

  private:
    const double min_density_;
    const double new_density_;
    const double new_pressure_;
    const EquationOfState& eos_;
  };
}

int main(void)
{
  // Set up the boundary type for the points
  SquareBox outer(Vector2D(-1,-1), Vector2D(1,1));

  // Set up the initial grid points
  int npointsx=50;
  int npointsy=50;
  vector<Vector2D> InitPoints=
    cartesian_mesh(npointsx,npointsy,
		   outer.getBoundary().first,
		   outer.getBoundary().second);

  // Set up the tessellation
  VoronoiMesh tess;

  // Set up the Riemann solver
  Hllc rs;

  // Set the hydro boundary conditions
  RigidWallHydro hbc(rs);

  // Set up the equation of state
  double gamma=5./3.;
  IdealGas eos(gamma);

  // Set up the point motion scheme
  Lagrangian l_motion;
  RoundCells pointmotion(l_motion,hbc);

  // Set up the interpolation
  LinearGaussConsistent interpolation(eos,outer,hbc);

  // Set up the initial Hydro
  Uniform2D density(1);
  Circle hot_spot(Vector2D(0,0),0.3);
  Uniform2D low_pressure(1);
  Uniform2D high_pressure(100);
  Piecewise pressure(hot_spot, high_pressure, low_pressure);
  Uniform2D xvelocity(0);
  Uniform2D yvelocity(0);

  // Set up the external source term
  ZeroForce force;

  // Set up the simulation
  hdsim sim(InitPoints,tess,interpolation,density,pressure,xvelocity,
	    yvelocity,eos,rs,pointmotion,force,outer,hbc);

  // Choose the Courant number
  sim.SetCfl(0.7);

  // How long shall we run the simulation?
  const double tend=0.05;
  sim.SetEndTime(tend);

  // Main simulation loop
  SafeTimeTermination term_cond(tend,1e6);
  ProgressReport diag1(25);
  ConsecutiveSnapshots diag2(0.01);
  MultipleDiagnostics diag;
  diag.diag_list.push_back(&diag1);
  diag.diag_list.push_back(&diag2);
  CellHeater cell_heater(0.5,1,10,eos);
  simulation2d::main_loop(sim,
			  term_cond,
			  &hdsim::TimeAdvance2Mid,
			  &diag,
			  &cell_heater);

  // Done running the simulation, output the data
  write_snapshot_to_hdf5(sim,"final.h5");

  // We are done!!
  cout<<"Finished running the simulation"<<endl;

  return 0;
}
