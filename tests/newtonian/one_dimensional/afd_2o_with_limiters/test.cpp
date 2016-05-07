#include <iostream>
#include <cmath>
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/plm1d.hpp"
#include "source/newtonian/one_dimensional/eos_consistent1d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/rigid_wall_1d.hpp"
#include "source/newtonian/one_dimensional/zero_force_1d.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/test_1d/simple_waves_ideal_gas.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_1d/main_loop_1d.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"

using namespace std;
using namespace simulation1d;
using namespace diagnostics1d;
using namespace interpolations1d;

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
    edges_(linspace(0,1,500)),
    eos_(read_number("adiabatic_index.txt")),
    naive_(),
    interpm_(naive_, eos_),
    rs_(),
    vm_(),
    bc_(),
    force_(),
    sim_(pg_,
	 edges_,
	 interpm_,
	 Uniform(read_number("ambient_density.txt")),
	 GaussianPert(read_number("ambient_pressure.txt"),
		      1e-3,
		      0.01,
		      0.5),
	 Uniform(read_number("ambient_velocity.txt")),
	 Uniform(0),
	 eos_,
	 rs_, 
	 vm_, 
	 bc_,
	 force_) {}

  hdsim1D& getSim(void)
  {
    return sim_;
  }
  
private:
  const SlabSymmetry1D pg_;
  vector<double> edges_;
  IdealGas eos_;
  PLM1D naive_;
  EOSConsistent interpm_;
  Hllc rs_;
  Eulerian1D vm_;
  RigidWall1D bc_;
  const ZeroForce1D force_;
  hdsim1D sim_;
};

int main(void)
{
  SimData sim_data;
  hdsim1D& sim = sim_data.getSim();
  
  write_snapshot_to_hdf5(sim,"initial.h5");

  main_loop(sim, 0.3, 1e6, 2,
	    "time.txt");

  write_snapshot_to_hdf5(sim,"final.h5");

  return 0;
}
