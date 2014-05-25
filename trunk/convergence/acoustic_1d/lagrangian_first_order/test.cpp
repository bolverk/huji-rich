#include <iostream>
#include <fstream>
#include <cmath>
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/one_dimensional/plm1d.hpp"
#include "source/newtonian/one_dimensional/eos_consistent1d.hpp"
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/lagrangian1d.hpp"
#include "source/newtonian/one_dimensional/periodic_1d.hpp"
#include "source/newtonian/one_dimensional/diagnostics_1d.hpp"
#include "source/misc/int2str.hpp"
#include "source/misc/utils.hpp"
#include "source/convergence/external_data.hpp"
#include "source/convergence/choose_between.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/one_dimensional/zero_force_1d.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"
#include "source/newtonian/test_1d/main_loop_1d.hpp"
#include "source/newtonian/test_1d/acoustic.hpp"

using namespace std;
using namespace interpolations1d;
using namespace diagnostics1d;
using namespace simulation1d;

class SimData
{
public:

  /*! \brief Class constructor
    \param exd Data from external file
   */
  SimData(int res):
    vertices_(linspace(0,1,res)),
    eos_(5./3.),
    pcm_(),
    init_cond_(1,3./5.,eos_,1e-6,1),
    rs_(),
    lagrangian_(true),
    bc_(),
    force_(),
    sim_(vertices_,
	 pcm_,
	 init_cond_.getProfile("density"),
	 init_cond_.getProfile("pressure"),
	 init_cond_.getProfile("xvelocity"),
	 init_cond_.getProfile("yvelocity"),
	 eos_,
	 rs_,
	 lagrangian_,
	 bc_,
	 force_)
  {}

  /*! \brief Returns a reference to the simulation
   */
  hdsim1D& getSim(void)
  {
    return sim_;
  }

private:
  const vector<double> vertices_;
  const IdealGas eos_;
  const PCM1D pcm_;
  const AcousticInitCond init_cond_;
  const Hllc rs_;
  const Lagrangian1D lagrangian_;
  const Periodic1D bc_;
  const ZeroForce1D force_;
  hdsim1D sim_;
};

int main(void)
{
  SimData sim_data(read_int("resolution.txt"));
  hdsim1D& sim = sim_data.getSim();
	
  main_loop(sim, 1, 1e6, 1,
	    "time.txt");

  write_snapshot_to_hdf5(sim,"final.h5");
  
  return 0;
}
