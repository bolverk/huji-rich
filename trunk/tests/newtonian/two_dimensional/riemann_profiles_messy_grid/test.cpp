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
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/test_2d/square_grid.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"

using namespace std;
using namespace simulation2d;

class SimData
{
public:

  SimData(void):
    width_(1),
    init_points_(square_grid(width_,30)),
    outer_(0,width_,width_,0),
    tess_(),
    interp_method_(),
    density_(1),
    pressure_(0,0.5,0,1,2,1),
    xvelocity_(0),
    yvelocity_(0),
    eos_(5./3.),
    pm_naive_(),
    rs_(),
    hbc_(rs_),
    point_motion_(pm_naive_,hbc_),
    force_(),
    sim_(init_points_,
	 tess_,
	 interp_method_,
	 density_,
	 pressure_,
	 xvelocity_,
	 yvelocity_,
	 eos_,
	 rs_,
	 point_motion_,
	 force_,
	 outer_,
	 hbc_) {}

  hdsim& getSim(void)
  {
    return sim_;
  }

private:
  const double width_;
  const vector<Vector2D> init_points_;
  const SquareBox outer_;
  VoronoiMesh tess_;
  PCM2D interp_method_;
  const Uniform2D density_;
  const Step2D pressure_;
  const Uniform2D xvelocity_;
  const Uniform2D yvelocity_;
  const IdealGas eos_;
  Lagrangian pm_naive_;
  const Hllc rs_;
  const RigidWallHydro hbc_;
  RoundCells point_motion_;
  ZeroForce force_;
  hdsim sim_;
};

namespace {
  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.225, 1e6);
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

  my_main_loop(sim);

  write_snapshot_to_hdf5(sim,"final.h5");

  return 0;
}

