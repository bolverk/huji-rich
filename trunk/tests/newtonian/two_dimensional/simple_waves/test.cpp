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
#include "source/newtonian/two_dimensional/interpolations/pcm2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/step2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/test_2d/square_grid.hpp"
#include "source/newtonian/test_2d/profile_1d.hpp"
#include "source/newtonian/test_1d/collela.hpp"
#include "source/newtonian/test_1d/simple_waves_ideal_gas.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"

using namespace std;
using namespace simulation2d;

namespace {
class InitProfiles
{
public:

  InitProfiles(void):
    width_(1),
    g_(read_number("adiabatic_index.txt")),
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

  SimData(InitProfiles init_prof = InitProfiles()):
    outer_(0,init_prof.getWidth(),
	   init_prof.getWidth(),0),
    tess_(),
    interpm_(),
    eos_(init_prof.getAdiabaticIndex()),
    bpm_(),
    rs_(),
    hbc_(rs_),
    point_motion_(bpm_,hbc_),
    force_(),
    sim_(square_grid(init_prof.getWidth(),30),
	 &tess_,
	 &interpm_,
	 init_prof.getDensity(),
	 init_prof.getPressure(),
	 init_prof.getXVelocity(),
	 init_prof.getYVelocity(),
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
  const SquareBox outer_;
  VoronoiMesh tess_;
  PCM2D interpm_;
  const IdealGas eos_;
  Lagrangian bpm_;
  const Hllc rs_;
  const RigidWallHydro hbc_;
  RoundCells point_motion_;
  ZeroForce force_;
  hdsim sim_;
};

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.02, 1e6);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      1,
	      &diag);
  }
}

int main(void)
{
  SimData sim_data;
  hdsim& sim = sim_data.getSim();

  write_snapshot_to_hdf5(sim, "initial.h5");

  my_main_loop(sim);

  write_snapshot_to_hdf5(sim, "final.h5");

  return 0;
}

