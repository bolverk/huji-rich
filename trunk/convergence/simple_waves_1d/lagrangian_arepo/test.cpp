#include <iostream>
#include <fstream>
#include <cmath>
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/rigid_wall_1d.hpp"
#include "source/newtonian/one_dimensional/diagnostics_1d.hpp"
#include "source/misc/utils.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/one_dimensional/plm1d.hpp"
#include "source/newtonian/one_dimensional/lagrangian1d.hpp"
#include "source/newtonian/one_dimensional/zero_force_1d.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"
#include "source/newtonian/test_1d/main_loop_1d.hpp"
#include "source/newtonian/test_1d/collela.hpp"
#include "source/newtonian/test_1d/simple_waves_ideal_gas.hpp"
#include "source/newtonian/one_dimensional/arepo_interp.hpp"

using namespace std;
using namespace diagnostics1d;
using namespace simulation1d;

/*! \brief Initial conditions giving rise to simple waves
 */
class SimpleWaveInitCond
{
public:

  SimpleWaveInitCond
  (double ref, double a,
   double l, double offset,
   double s, double g):
    density_(ref,a,l,offset),
    pressure_(density_,s,g),
    eos_(g),
    sound_speed_(pressure_,
		 density_,
		 eos_),
    c_edge_(sound_speed_.EvalAt(offset+l+1)),
    rv_edge_(calc_riemann_invariant(0,c_edge_,g,0)),
    xvelocity_(rv_edge_,0,sound_speed_,g),
    yvelocity_(0) {}

  IdealGas const& getEOS(void) const
  {
    return eos_;
  }

  SpatialDistribution1D const& getProfile(string pname) const
  {
    if("density"==pname)
      return density_;
    else if("pressure"==pname)
      return pressure_;
    else if("xvelocity"==pname)
      return xvelocity_;
    else if("yvelocity"==pname)
      return yvelocity_;
    else
      throw "Unknown variable name "+pname;
  }

private:
  
  const Collela density_;
  const ConstEntropy pressure_;
  const IdealGas eos_;
  const SoundSpeedDist sound_speed_;
  const double c_edge_;
  const double rv_edge_;
  const ConstRiemannInv xvelocity_;
  const Uniform yvelocity_;
};

class SimData
{
public:
  
  SimData(int res):
    edges_(linspace(0,1,res)),
    init_cond_(1,10,0.3,0.5,1,5./3.),
    interp_(init_cond_.getEOS()),
    rs_(),
    lagrangian_(true),
    bc_(),
    force_(),
    sim_(edges_,
	 interp_,
	 init_cond_.getProfile("density"),
	 init_cond_.getProfile("pressure"),
	 init_cond_.getProfile("xvelocity"),
	 init_cond_.getProfile("yvelocity"),
	 init_cond_.getEOS(),
	 rs_,
	 lagrangian_,
	 bc_,
	 force_) {}

  hdsim1D& getSim(void)
  {
    return sim_;
  }
  
private:
  const vector<double> edges_;
  const SimpleWaveInitCond init_cond_;
  const ArepoInterp interp_;
  const Hllc rs_;
  const Lagrangian1D lagrangian_;
  const RigidWall1D bc_;
  const ZeroForce1D force_;
  hdsim1D sim_;
};

int main(void)
{
  SimData sim_data(read_int("resolution.txt"));
  hdsim1D& sim = sim_data.getSim();
  
  write_snapshot_to_hdf5(sim,"initial.h5");

  main_loop(sim, 0.02, 1e6, 1,
	    "time.txt");

  write_snapshot_to_hdf5(sim,"final.h5");

  return 0;
}
