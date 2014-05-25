#include <iostream>
#include "source/newtonian/test_2d/square_grid.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/interpolations/pcm2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/Circle2D.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylinderical_geometry.hpp"
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
    offset_(1e-3,0),
    mesh_(offset_grid(square_grid(1,50),offset_)),
    outer_(offset_.x,offset_.x+1,
	   offset_.y+1,offset_.y),
    tess_(),
    interpm_(),
    density_(1),
    pressure_(0+offset_.x,
	      0.5+offset_.y,
	      0.05,1e5,1),
    xvelocity_(0),
    yvelocity_(0),
    eos_(5./3.),
    rs_(),
    hbc_(rs_),
    raw_point_motion_(),
    point_motion_(raw_point_motion_),
    force_(Vector2D(0,0),
	   Vector2D(0,1)),
    sim_(mesh_,
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
  const Vector2D offset_;
  const vector<Vector2D> mesh_;
  const SquareBox outer_;
  VoronoiMesh tess_;
  PCM2D interpm_;
  const Uniform2D density_;
  const Circle2D pressure_;
  const Uniform2D xvelocity_;
  const Uniform2D yvelocity_;
  const IdealGas eos_;
  const Hllc rs_;
  const RigidWallHydro hbc_;
  Lagrangian raw_point_motion_;
  RoundCells point_motion_;
  CylindericalGeometry force_;
  hdsim sim_;
};

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.01,1e6);
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

  write_snapshot_to_hdf5(sim, "final.h5");

  return 0;
}
