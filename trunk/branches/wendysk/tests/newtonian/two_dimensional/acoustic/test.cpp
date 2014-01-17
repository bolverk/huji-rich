#include <iostream>
#include "source/newtonian/test_1d/acoustic.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/PeriodicBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/PeriodicHydro.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/interpolations/linear_gauss.hpp"
#include "source/newtonian/two_dimensional/interpolations/eos_consistent.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/test_2d/profile_1d.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/test_2d/square_grid.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"

using namespace std;
using namespace simulation2d;

namespace {

class SimData
{
public:

  SimData(void):
    width_(read_number("width.txt")),
    init_points_(square_grid(width_,30)),
    outer_(0,width_,width_,0),
    tess_(),
    interpm_naive_(outer_,&hbc_,true,false),
    eos_(read_number("adiabatic_index.txt")),
    interpm_(interpm_naive_,eos_),
    init_cond_(read_number("ambient_density.txt"),
	       read_number("ambient_pressure.txt"),
	       eos_,
	       read_number("amplitude.txt"),
	       width_),
    density_(init_cond_.getProfile("density")),
    pressure_(init_cond_.getProfile("pressure")),
    xvelocity_(init_cond_.getProfile("xvelocity")),
    yvelocity_(init_cond_.getProfile("yvelocity")),
    bpm_(),
    point_motion_(bpm_),
    rs_(),
    hbc_(rs_),
    force_(),
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

  double width_;
  vector<Vector2D> init_points_;
  PeriodicBox outer_;
  VoronoiMesh tess_;
  LinearGauss interpm_naive_;
  IdealGas eos_;
  EOSConsistent interpm_;
  AcousticInitCond init_cond_;
  Profile1D density_;
  Profile1D pressure_;
  Profile1D xvelocity_;
  Profile1D yvelocity_;
  Lagrangian bpm_;
  RoundCells point_motion_;
  Hllc rs_;
  PeriodicHydro hbc_;
  ZeroForce force_;
  hdsim sim_;
};
}

int main(void)
{
  SimData sim_data;
  hdsim& sim = sim_data.getSim();

  {
    SafeTimeTermination term_cond(1,1e6);
    WriteTime diag("time.txt");
    main_loop(sim, 
	      term_cond,
	      2,
	      &diag);
  }

  write_snapshot_to_hdf5(sim, "final.h5");

  return 0;
}
