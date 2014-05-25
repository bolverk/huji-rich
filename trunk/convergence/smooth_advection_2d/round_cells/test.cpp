#include <iostream>
#include <cmath>
#include "source/newtonian/one_dimensional/diagnostics_1d.hpp"
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/two_dimensional/PeriodicBox.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/pcm2d.hpp"
#include "source/newtonian/two_dimensional/uniform2d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/eulerian.hpp"
#include "source/newtonian/two_dimensional/lagrangian.hpp"
#include "source/newtonian/two_dimensional/round_cells.hpp"
#include "source/newtonian/two_dimensional/PeriodicHydro.hpp"
#include "source/newtonian/two_dimensional/zero_force.hpp"
#include "source/newtonian/test_2d/square_grid.hpp"
#include "source/newtonian/test_2d/profile_1d.hpp"
#include "source/misc/simple_io.hpp"

using namespace std;

class SineWaveCellAverage: public SpatialDistribution1D
{
public:

  SineWaveCellAverage(double amplitude,
		      double wavelength,
		      double phase,
		      double offset,
		      double cell_width):
    amp_(amplitude),
    k_(2*M_PI/wavelength),
    ph_(phase),
    offset_(offset),
    dx_(cell_width) {}

  double EvalAt(double x) const
  {
    return offset_ + 
      (amp_/(dx_*k_))*
      (cos(k_*(x-dx_/2))-
       cos(k_*(x+dx_/2)));
  }

private:

  const double amp_;
  const double k_;
  const double ph_;
  const double offset_;
  const double dx_;
};

class SineWave: public SpatialDistribution1D
{
public:

  SineWave(double amplitue,
	   double wavelength,
	   double phase,
	   double offset):
    amp_(amplitue),
    k_(2*M_PI/wavelength),
    phase_(phase),
    offset_(offset) {}

  double EvalAt(double x) const
  {
    return amp_*sin(k_*x+phase_)+offset_;
  }

private:

  const double amp_;
  const double k_;
  const double phase_;
  const double offset_;
};

class SimData
{
public:

  SimData(double width, int resolution):
    init_points_(square_grid(width,resolution)),
    outer_(0,width,width,0),
    eos_(5./3.),
    density_1d_(1,width,0,2),
    density_(density_1d_),
    pressure_(1),
    xvelocity_(1),
    yvelocity_(0),
    point_motion_(naive_pm_),
    hbc_(rs_),
    sim_(init_points_,
	 &tess_,
	 &interpm_,
	 density_,
	 pressure_,
	 xvelocity_,
	 yvelocity_,
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
  const vector<Vector2D> init_points_;
  const PeriodicBox outer_;
  VoronoiMesh tess_;
  const IdealGas eos_;
  const SineWave density_1d_;
  const Profile1D density_;
  const Uniform2D pressure_;
  const Uniform2D xvelocity_;
  const Uniform2D yvelocity_;
  PCM2D interpm_;
  const Hllc rs_;
  Lagrangian naive_pm_;
  RoundCells point_motion_;
  const PeriodicHydro hbc_;
  ZeroForce force_;
  hdsim sim_;
};

void main_loop(hdsim& sim)
{
  const double tf = 1;
  const int max_cycle = (int)1e6;

  while(sim.GetTime()<tf){
    sim.TimeAdvance();

    write_number(sim.GetTime(),"time.txt");

    if(sim.GetCycle()>max_cycle){
      cout << "Max number of cycles exceeded in main loop" << endl;
      throw;
    }
  }
}

void write_output(hdsim const& sim)
{
  write_x_plot(sim,"sim_results.txt");
}

int main(void)
{
  SimData sim_data(1,read_int("resolution.txt"));
  hdsim& sim = sim_data.getSim();

  main_loop(sim);

  write_output(sim);

  return 0;
}
