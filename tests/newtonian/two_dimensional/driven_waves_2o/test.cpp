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
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/ConservativeForce.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/PeriodicBox.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/modular_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/ghost_point_generators/PeriodicGhostGenerator.hpp"
#include "source/newtonian/two_dimensional/interpolations/LinearGaussImproved.hpp"
#include "source/newtonian/two_dimensional/periodic_edge_velocities.hpp"

using namespace std;
using namespace simulation2d;

namespace {
  vector<ComputationalCell> calc_init_cond
  (const Tessellation& tess)
  {
    vector<ComputationalCell> res
      (static_cast<size_t>(tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i){
      res[i].density = 1;
      res[i].pressure = 1;
      res[i].velocity = Vector2D(0.1,0);
    }
    return res;
  }

  class PeriodicDriver: public Acceleration
  {
  public:

    PeriodicDriver(double wavelenth,
		   double amplitude):
      k_(2*M_PI/wavelenth),
      amp_(amplitude) {}

    Vector2D operator()
    (const Tessellation& tess,
     const vector<ComputationalCell>& /*cells*/,
     const vector<Extensive>& /*fluxes*/,
     double /*t*/,
     int point,
	 TracerStickerNames const& /*ts*/) const
    {
      return Vector2D(-amp_*cos(k_*tess.GetCellCM(point).x),0);
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
      init_points_(cartesian_mesh(30,30,Vector2D(0,0),
				  Vector2D(width_,width_))),
      outer_(0,width_,width_,0),
      tess_(init_points_, outer_),
      density_(1),
      pressure_(1),
      xvelocity_(0.1),
      yvelocity_(0),
      eos_(5./3.),
      pm_naive_(),
      rs_(),
      point_motion_(pm_naive_,eos_),
      evc_(),
      acc_(1,0.001),
      force_(acc_),
      tsf_(0.3),
      gpg_(),
      sr_(eos_,gpg_),
      fc_(sr_,rs_),
      eu_(),
      sim_(tess_,
	   outer_,
	   pg_,
	   calc_init_cond(tess_),
	   eos_,
	   point_motion_,
	   evc_,
	   force_,
	   tsf_,
	   fc_,
	   eu_,
	   cu_) {}

    hdsim& getSim(void)
    {
      return sim_;
    }

  private:
    const SlabSymmetry pg_;
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
    RoundCells point_motion_;
    const PeriodicEdgeVelocities evc_;
    PeriodicDriver acc_;
    ConservativeForce force_;
    const SimpleCFL tsf_;
    const PeriodicGhostGenerator gpg_;
    const LinearGaussImproved sr_;
    const ModularFluxCalculator fc_;
    const SimpleExtensiveUpdater eu_;
    const SimpleCellUpdater cu_;
    hdsim sim_;
  };
  
  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(1, 1e6);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      &hdsim::TimeAdvance2Heun,
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

