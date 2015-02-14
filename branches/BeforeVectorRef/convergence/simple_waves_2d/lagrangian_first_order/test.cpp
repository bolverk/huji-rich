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
#include "source/newtonian/two_dimensional/uniform2d.hpp"
#include "source/newtonian/two_dimensional/step2d.hpp"
#include "source/newtonian/two_dimensional/eulerian.hpp"
#include "source/newtonian/two_dimensional/lagrangian.hpp"
#include "source/newtonian/two_dimensional/SquareBox.hpp"
#include "source/newtonian/two_dimensional/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/round_cells.hpp"
#include "source/newtonian/two_dimensional/zero_force.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/test_2d/square_grid.hpp"
#include "source/newtonian/test_2d/profile_1d.hpp"
#include "source/newtonian/test_1d/collela.hpp"
#include "source/newtonian/test_1d/simple_waves_ideal_gas.hpp"
#include "source/misc/simple_io.hpp"

using namespace std;

class InitProfiles
{
public:

  InitProfiles(void):
    width_(1),
    g_(5./3.),
    density_1d_(1,10,0.3,0.5),
    swigic_(density_1d_,1,g_,width_),
    density_(density_1d_),
    pressure_(swigic_.getProfile("pressure")),
    xvelocity_(swigic_.getProfile("xvelocity")),
    yvelocity_(0) {}

  double getWidth(void)
  {
    return width_;
  }

  IdealGas const& getEOS(void) const
  {
    return swigic_.getEOS();
  }

  double getAdiabaticIndex(void) const
  {
    return g_;
  }

  Profile1D const& getDensity(void) const
  {
    return density_;
  }

  Profile1D const& getPressure(void) const
  {
    return pressure_;
  }

  Profile1D const& getXVelocity(void) const
  {
    return xvelocity_;
  }

  Uniform2D const& getYVelocity(void) const
  {
    return yvelocity_;
  }

private:
  const double width_;
  const double g_;
  const Collela density_1d_;
  const SimpleWaveIdealGasInitCond swigic_;
  const Profile1D density_;
  const Profile1D pressure_;
  const Profile1D xvelocity_;
  const Uniform2D yvelocity_;
};

class SimData 
{
public:

  SimData(int res, InitProfiles init_prof = InitProfiles()):
    outer_(0,init_prof.getWidth(),
	   init_prof.getWidth(),0),
    eos_(init_prof.getAdiabaticIndex()),
    point_motion_(bpm_),
    hbc_(rs_),
    sim_(square_grid(init_prof.getWidth(),res),
	 &tess_,
	 &interpm_,
	 init_prof.getDensity(),
	 init_prof.getPressure(),
	 init_prof.getXVelocity(),
	 init_prof.getYVelocity(),
	 eos_,
	 &rs_,
	 &point_motion_,
	 &force_,
	 &outer_,
	 &hbc_) {}

  hdsim& getSim(void)
  {
    return sim_;
  }

private:
  const SquareBox outer_;
  VoronoiMesh tess_;
  PCM2D interpm_;
  const IdealGas eos_;
  Lagrangian bpm_;
  RoundCells point_motion_;
  const Hllc rs_;
  const RigidWallHydro hbc_;
  ZeroForce force_;
  hdsim sim_;
};

void main_loop(hdsim& sim)
{
  const int max_iter = 5e6;
  const double tf = 0.02;
  while(tf>sim.GetTime()){

    sim.TimeAdvance();
    
    write_number(sim.GetTime(),"time.txt");

    if(sim.GetCycle()>max_iter)
      throw "Maximum number of iterations exceeded in main loop";
  }
}

void write_output(hdsim const& sim)
{
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
  write_x_plot(sim,"x_prof_final.txt");
}

int main(void)
{
  SimData sim_data(read_int("resolution.txt"));
  hdsim& sim = sim_data.getSim();

  write_x_plot(sim,"x_prof_initial.txt");

  main_loop(sim);

  write_output(sim);

  return 0;
}

