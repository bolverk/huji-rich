#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "source/tessellation/geometry.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/pcm2d.hpp"
#include "source/newtonian/two_dimensional/linear_gauss.hpp"
#include "source/newtonian/two_dimensional/uniform2d.hpp"
#include "source/newtonian/two_dimensional/step2d.hpp"
#include "source/newtonian/two_dimensional/lagrangian.hpp"
#include "source/newtonian/two_dimensional/round_cells.hpp"
#include "source/newtonian/two_dimensional/zero_force.hpp"
#include "source/newtonian/two_dimensional/SquareBox.hpp"
#include "source/newtonian/two_dimensional/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/test_2d/square_grid.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/PeriodicBox.hpp"
#include "source/newtonian/two_dimensional/PeriodicHydro.hpp"
#include "source/newtonian/two_dimensional/eos_consistent.hpp"

using namespace std;

class PeriodicDriver: public SourceTerm
{
public:

  PeriodicDriver(double wavelenth,
		 double amplitude,
		 double phase_velocity):
    k_(2*M_PI/wavelenth),
    amp_(amplitude),
    v_(phase_velocity) {}

  Conserved Calculate
  (Tessellation const* tess,
   vector<Primitive> const& cells,
   int point,
   vector<Conserved> const& /*fluxes*/,
   vector<Vector2D> const& /*point_velocity*/,
   HydroBoundaryConditions const* /*hbc*/,
   double t,
   double /*dt*/)
  {
    const double volume = tess->GetVolume(point);
    const double x = tess->GetMeshPoint(point).x;
    const double density = cells[point].Density;
    const double xvelocity = cells[point].Velocity.x;
    const double acceleration = amp_*sin(k_*x)*sin(k_*v_*t);
    const double xmom = density*acceleration;
    const double enr = density*acceleration*xvelocity;
    return -volume*Conserved
      (0,Vector2D(xmom,0),enr);
  }

private:
  const double k_;
  const double amp_;
  const double v_;
};

class SimData
{
public:

  SimData(int res):
    width_(1),
    init_points_(square_grid(width_,res)),
    outer_(0,width_,width_,0),
    density_(read_number("mean_density.txt")),
    pressure_(read_number("mean_pressure.txt")),
    xvelocity_(0),
    yvelocity_(0),
    eos_(read_number("adiabatic_index.txt")),
    point_motion_(pm_naive_),
    hbc_(rs_),
    interp_naive_(outer_,&hbc_,true,false),
    interpm_(interp_naive_,eos_),
    force_(read_number("wavelength.txt"),
	   read_number("amplitude.txt"),
	   read_number("phase_velocity.txt")),
    sim_(init_points_,
	 &tess_,
	 &interpm_,
	 density_,
	 pressure_,
	 xvelocity_,
	 yvelocity_,
	 eos_,
	 rs_,
	 &point_motion_,
	 &force_,
	 &outer_,
	 &hbc_) {}

  hdsim& getSim(void)
  {
    return sim_;
  }

private:
  const double width_;
  const vector<Vector2D> init_points_;
  const SquareBox outer_;
  VoronoiMesh tess_;
  //PCM2D interp_method_;
  const Uniform2D density_;
  const Uniform2D pressure_;
  const Uniform2D xvelocity_;
  const Uniform2D yvelocity_;
  const IdealGas eos_;
  Lagrangian pm_naive_;
  RoundCells point_motion_;
  const Hllc rs_;
  const RigidWallHydro hbc_;
  LinearGauss interp_naive_;
  EOSConsistent interpm_;
  PeriodicDriver force_;
  hdsim sim_;
};

void main_loop(hdsim& sim)
{
  const int max_iter = (int)1e6;
  const double tf = 1;
  
  while(sim.GetTime()<tf){
    sim.TimeAdvance2Mid();
    write_number(sim.GetTime(),"time.txt",14);

    if(sim.GetCycle()>max_iter)
      throw "Error in main_loop: max number of iterations exceeded";
  }
}

void write_output(hdsim const& sim)
{
  const int prec = 14;
  write_edges_and_neighbors(sim,"edges.txt");
  write_generating_points(sim,"pointpos.txt");
  write_cells_property(sim,"density",
		       "densities.txt");
  write_cells_property(sim,"pressure",
		       "pressures.txt");
  write_cells_property(sim,"velocity x",
		       "xvelocities.txt");
  write_cells_property(sim,"velocity y",
		       "yvelocities.txt");
  write_x_plot(sim,"prof1d.txt",prec);
}

int main(void)
{
  SimData sim_data(read_int("resolution.txt"));
  hdsim& sim = sim_data.getSim();

  main_loop(sim);
  
  write_output(sim);

  return 0;
}

