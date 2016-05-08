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
#include "source/newtonian/test_1d/collela.hpp"
#include "source/newtonian/test_1d/simple_waves_ideal_gas.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_1d/main_loop_1d.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"
#include "source/newtonian/one_dimensional/plm1d.hpp"
#include "source/newtonian/one_dimensional/eos_consistent1d.hpp"

using namespace std;
using namespace simulation1d;
using namespace diagnostics1d;
using namespace interpolations1d;

/*! \brief Spatial profiles giving rise only to a simple wave, with collela distribution for the density
 */
class CollelaSimpleWave
{
public:

  CollelaSimpleWave(double ref, double a,
		    double l, double offset,
		    double s, double g):
    density_(ref,a,l,offset),
    init_cond_(density_,s,g) {}

  IdealGas const& getEOS(void) const
  {
    return init_cond_.getEOS();
  }

  SpatialDistribution1D const& getProfile(string const& pname) const
  {
    return init_cond_.getProfile(pname);
  }

private:

  Collela density_;
  SimpleWaveIdealGasInitCond init_cond_;

};

class SimData
{
public:
  
  SimData(void):
    edges_(linspace(0,1,100)),
    init_cond_(1,10,0.3,0.5,1,5./3.),
    naive_(),
    interpm_(naive_, init_cond_.getEOS()),
    rs_(),
    vm_(),
    bc_(),
    force_(),
    sim_
    (pg_,
     edges_,
     interpm_,
     init_cond_.getProfile("density"),
     init_cond_.getProfile("pressure"),
     init_cond_.getProfile("xvelocity"),
     init_cond_.getProfile("yvelocity"),
     init_cond_.getEOS(),
     rs_, vm_, bc_,
     force_) {}

  hdsim1D& getSim(void)
  {
    return sim_;
  }
  
private:
  const SlabSymmetry1D pg_;
  const vector<double> edges_;
  CollelaSimpleWave init_cond_;
  PLM1D naive_;
  EOSConsistent interpm_;
  const Hllc rs_;
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

  main_loop(sim, 0.02, 1e6, 1,
	    "time.txt");

  write_snapshot_to_hdf5(sim,"final.h5");

  return 0;
}
