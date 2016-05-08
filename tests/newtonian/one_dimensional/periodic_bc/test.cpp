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

using namespace std;
using namespace simulation1d;
using namespace diagnostics1d;

namespace {
class SimData
{
public:

  SimData(void):
    vertices_(linspace(0,1,1000)),
    interpm_(),
    eos_(5./3.),
    init_cond_(1,3./5.,eos_,1e-6,1),
    rs_(),
    vm_(),
    bc_(),
    force_(),
    sim_(pg_,
	 vertices_,
	 interpm_,
	 init_cond_.getProfile("density"),
	 init_cond_.getProfile("pressure"),
	 init_cond_.getProfile("xvelocity"),
	 init_cond_.getProfile("yvelocity"),
	 eos_,
	 rs_,
	 vm_,
	 bc_,
	 force_)
  {}

  hdsim1D& getSim(void)
  {
    return sim_;
  }

private:
  const SlabSymmetry1D pg_;
  const vector<double> vertices_;
  PCM1D interpm_;
  const IdealGas eos_;
  const AcousticInitCond init_cond_;
  const Hllc rs_;
  const Eulerian1D vm_;
  const Periodic1D bc_;
  const ZeroForce1D force_;
  hdsim1D sim_;
};
}

int main(void)
{
  SimData sim_data;
  hdsim1D& sim = sim_data.getSim();
			
  main_loop(sim, 1, 1e6, 1, "time.txt");

  write_snapshot_to_hdf5(sim,"final.h5");
  
  return 0;
}
