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
#include "source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/step2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/ConservativeForce.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/test_2d/square_grid.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/PeriodicBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/PeriodicHydro.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"

using namespace std;
using namespace simulation2d;

class PeriodicDriver: public Acceleration
{
public:

  PeriodicDriver(double wavelenth,
		 double amplitude):
    k_(2*M_PI/wavelenth),
    amp_(amplitude) {}

  Vector2D Calculate
  (Tessellation const* tess,
   vector<Primitive> const& /*cells*/,
   int point,
   vector<Conserved> const& /*fluxes*/,
   vector<Vector2D> const& /*point_velocity*/,
   HydroBoundaryConditions const* /*hbc*/,
   double /*t*/,
   double /*dt*/)
  {
    return Vector2D(-amp_*cos(k_*tess->GetCellCM(point).x),0);
  }

private:
  const double k_;
  const double amp_;
};

class SimData
{
public:

  SimData(void):
    width_(1),
    init_points_(square_grid(width_,30)),
    outer_(0,width_,width_,0),
    tess_(),
    density_(1),
    pressure_(1),
    xvelocity_(0.1),
    yvelocity_(0),
    eos_(5./3.),
    pm_naive_(),
    rs_(),
    hbc_(rs_),
    point_motion_(pm_naive_,hbc_),
    interpm_(eos_,outer_,&hbc_,true,false),
    acc_(1,0.001),
    force_(&acc_),
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
  const PeriodicBox outer_;
  VoronoiMesh tess_;
  const Uniform2D density_;
  const Uniform2D pressure_;
  const Uniform2D xvelocity_;
  const Uniform2D yvelocity_;
  const IdealGas eos_;
  Lagrangian pm_naive_;
  const Hllc rs_;
  const PeriodicHydro hbc_;
  RoundCells point_motion_;
  LinearGaussConsistent interpm_;
  PeriodicDriver acc_;
  ConservativeForce force_;
  hdsim sim_;
};

namespace {
  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(1, 1e6);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      2,
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

