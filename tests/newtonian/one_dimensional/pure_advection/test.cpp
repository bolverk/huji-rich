#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/assign/list_of.hpp>
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/one_dimensional/plm1d.hpp"
#include "source/newtonian/one_dimensional/eos_consistent1d.hpp"
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/periodic_1d.hpp"
#include "source/newtonian/one_dimensional/zero_force_1d.hpp"
#include "source/misc/int2str.hpp"
#include "source/misc/utils.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"
#include "source/newtonian/test_1d/main_loop_1d.hpp"

using namespace std;
using namespace interpolations1d;
using namespace simulation1d;
using namespace diagnostics1d;

namespace {
SpatialReconstruction1D const& choose_between
(SpatialReconstruction1D const& opt1,
 SpatialReconstruction1D const& opt2,
 SpatialReconstruction1D const& opt3,
 string const& name1,
 string const& name2,
 string const& name3,
 string const& choice)
{
  if(name1==choice)
    return opt1;
  else if(name2==choice)
    return opt2;
  else if(name3==choice)
    return opt3;
  else
    throw "Unknown option "+choice;
}

class SimData
{
public:

  SimData(string const& interp_method_name):
    vertices_(linspace(0,1,100)),
    eos_(5./3.),
    pcm_(),
    plm_naive_(),
    plm_(plm_naive_,eos_),
    interpm_(choose_between(pcm_,plm_naive_,plm_,
			    "pcm","plm_naive","plm",
			    interp_method_name)),
    density_(1,0.3,2,0.7,1),
    pressure_(1),
    xvelocity_(1),
    yvelocity_(0),
    rs_(),
    vm_(),
    bc_(),
    force_(),
    sim_
    (pg_,
     vertices_,
     interpm_,
     density_,
     pressure_,
     xvelocity_,
     yvelocity_,
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
  const vector<double> vertices_;
  const IdealGas eos_;
  PCM1D pcm_;
  PLM1D plm_naive_;
  EOSConsistent plm_;
  SpatialReconstruction1D const& interpm_;
  const TwoSteps density_;
  const Uniform pressure_;
  const Uniform xvelocity_;
  const Uniform yvelocity_;
  const Hllc rs_;
  const Eulerian1D vm_;
  const Periodic1D bc_;
  const ZeroForce1D force_;
  hdsim1D sim_;
};
}

int main(void)
{
  const vector<string> interp_names = 
    boost::assign::list_of("pcm")("plm_naive")("plm");
  for(size_t i=0;i<interp_names.size();++i){
    SimData sim_data(interp_names[i]);

    main_loop(sim_data.getSim(), 1, 1e6, 1, "time.txt");

    write_snapshot_to_hdf5(sim_data.getSim(),
			   interp_names[i]+"_final.h5");
  }
  
  return 0;
}
