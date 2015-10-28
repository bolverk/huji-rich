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
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/test_2d/profile_1d.hpp"
#include "source/newtonian/test_1d/collela.hpp"
#include "source/newtonian/test_1d/simple_waves_ideal_gas.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"

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

  vector<ComputationalCell> calc_init_cond(const Tessellation& tess)
  {
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    const InitProfiles ip;
    for(size_t i=0;i<res.size();++i){
      const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
      res[i].density = ip.getDensity()(r);
      res[i].pressure = ip.getPressure()(r);
      res[i].velocity = Vector2D
	(ip.getXVelocity()(r),
	 ip.getYVelocity()(r));
    }
    return res;
  }

  class SimData 
  {
  public:

    SimData(InitProfiles init_prof = InitProfiles()):
      pg_(),
      outer_(0,init_prof.getWidth(),
	     init_prof.getWidth(),0),
      tess_(cartesian_mesh(30,30,
			  Vector2D(0,0),
			  Vector2D(init_prof.getWidth(),
				   init_prof.getWidth())),
	    outer_),
      eos_(init_prof.getAdiabaticIndex()),
      bpm_(),
      rs_(),
      point_motion_(bpm_,eos_),
      evc_(),
      force_(),
      tsf_(0.3),
      fc_(rs_),
      eu_(),
      cu_(),
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
    const SquareBox outer_;
    VoronoiMesh tess_;
    const IdealGas eos_;
    Lagrangian bpm_;
    const Hllc rs_;
    RoundCells point_motion_;
    const StationaryBox evc_;
    ZeroForce force_;
    const SimpleCFL tsf_;
    const SimpleFluxCalculator fc_;
    const SimpleExtensiveUpdater eu_;
    const SimpleCellUpdater cu_;
    hdsim sim_;
  };

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.02, 1e6);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      &hdsim::TimeAdvance,
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

