#include <iostream>
#include <fstream>
#include <cmath>
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/periodic_1d.hpp"
#include "source/misc/int2str.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/test_1d/sine_wave.hpp"
#include "source/newtonian/test_1d/acoustic.hpp"
#include "source/newtonian/one_dimensional/zero_force_1d.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"
#include "source/newtonian/test_1d/main_loop_1d.hpp"
#include "source/newtonian/one_dimensional/simple_cfl_1d.hpp"
#include "source/newtonian/one_dimensional/simple_extensive_updater_1d.hpp"
#include "source/newtonian/one_dimensional/simple_flux_calculator_1d.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/cfg/env.h"

using namespace std;
using namespace simulation1d;
using namespace diagnostics1d;

namespace {

  class Gaussian: public SpatialDistribution1D
  {
  public:

    Gaussian(const double mean,
	     const double std):
      mean_(mean), std_(std) {}

    double operator()(const double x) const
    {
      return exp(-pow((x-mean_)/std_,2));
    }

  private:

    const double mean_;
    const double std_;
  };

  class SimData
  {
  public:

    SimData(void):
      colour_func_(0.5,0.2),
      tracer_helper_(string("colour"), &colour_func_),
      vertices_(linspace(0,1,100)),
      interpm_(),
      eos_(5./3.),
      //    init_cond_(1,3./5.,eos_,1e-6,1),
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
	   (vertices_,
	    Uniform(1),
	    Uniform(1),
	    Uniform(1),
	    Uniform(0),
	    vector<pair<string, const SpatialDistribution1D*> >(1, tracer_helper_),
	    vector<pair<string, const BoolSpatialDistribution*> >()),
	   eos_,
	   vm_,
	   force_,
	   tsf_,
	   fc_,
	   eu_,
	   cu_)
    {}

    hdsim1D& getSim(void)
    {
      return sim_;
    }

  private:
    const Gaussian colour_func_;
    const pair<string, const SpatialDistribution1D*> tracer_helper_;
    const SlabSymmetry1D pg_;
    const vector<double> vertices_;
    PCM1D interpm_;
    const IdealGas eos_;
    const Hllc rs_;
    const Eulerian1D vm_;
    const Periodic1D bc_;
    const ZeroForce1D force_;
    const SimpleCFL1D tsf_;
    const SimpleFluxCalculator1D fc_;
    const SimpleExtensiveUpdater1D eu_;
    const SimpleCellUpdater1D cu_;
    hdsim1D sim_;
  };
}

int main(void)
{
  spdlog::cfg::load_env_levels();
  SimData sim_data;
  hdsim1D& sim = sim_data.getSim();

  write_snapshot_to_hdf5(sim,"initial.h5");
			
  main_loop(sim, 1, 1e6, 1, "time.txt");

  write_snapshot_to_hdf5(sim,"final.h5");
  
  return 0;
}
