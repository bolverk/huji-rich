#include <iostream>
#include <cmath>
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/rigid_wall_1d.hpp"
#include "source/newtonian/one_dimensional/zero_force_1d.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/test_1d/simple_waves_ideal_gas.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_1d/main_loop_1d.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"
#include "source/newtonian/one_dimensional/cell_updater_1d.hpp"
#include "source/newtonian/one_dimensional/simple_cfl_1d.hpp"
#include "source/newtonian/one_dimensional/simple_extensive_updater_1d.hpp"
#include "source/newtonian/one_dimensional/simple_flux_calculator_1d.hpp"

using namespace std;
using namespace simulation1d;
using namespace diagnostics1d;

namespace {
  class GaussianPert: public SpatialDistribution1D
  {
  public:

    GaussianPert(double ambient,
		 double amplitude,
		 double width,
		 double center):
      ambient_(ambient),
      amplitude_(amplitude),
      width_(width),
      center_(center) {}

    double operator()(double x) const
    {
      return ambient_ + amplitude_*exp(-pow((x-center_)/width_,2));
    }

  private:
    const double ambient_;
    const double amplitude_;
    const double width_;
    const double center_;
  };
}

class SimData
{
public:
  
  SimData(void):
    edges_(linspace(0,1,1000)),
    eos_(read_number("adiabatic_index.txt")),
    interpm_(),
    rs_(),
    vm_(),
    bc_(),
    force_(),
    tsf_(0.3),
    fc_(rs_,interpm_,bc_),
    eu_(),
    cu_(),
    sim_(pg_,
	 SimulationState1D
	 (edges_,
	  Uniform(read_number("ambient_density.txt")),
	  GaussianPert(read_number("ambient_pressure.txt"),
		       1e-3,
		       0.01,
		       0.5),
	  Uniform(read_number("ambient_velocity.txt")),
	  Uniform(0),
	  vector<pair<string, const SpatialDistribution1D*> >(),
	  vector<pair<string, const BoolSpatialDistribution*> >()),
	 eos_,
	 vm_, 
	 force_,
	 tsf_,
	 fc_,
	 eu_,
	 cu_) {}

  hdsim1D& getSim(void)
  {
    return sim_;
  }
  
private:
  const SlabSymmetry1D pg_;
  vector<double> edges_;
  IdealGas eos_;
  PCM1D interpm_;
  Hllc rs_;
  Eulerian1D vm_;
  RigidWall1D bc_;
  const ZeroForce1D force_;
  const SimpleCFL1D tsf_;
  const SimpleFluxCalculator1D fc_;
  const SimpleExtensiveUpdater1D eu_;
  const SimpleCellUpdater1D cu_;
  hdsim1D sim_;
};

int main(void)
{
  SimData sim_data;
  hdsim1D& sim = sim_data.getSim();
  
  write_snapshot_to_hdf5(sim,"initial.h5");

  main_loop(sim, 0.3, 1e6, 1,
	    "time.txt");

  write_snapshot_to_hdf5(sim,"final.h5");

  return 0;
}
