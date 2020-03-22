#include <iostream>
#include <fstream>
#include "source/newtonian/common/ersig.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/lagrangian1d.hpp"
#include "source/newtonian/one_dimensional/rigid_wall_1d.hpp"
#include "source/newtonian/one_dimensional/zero_force_1d.hpp"
#include "source/misc/int2str.hpp"
#include "source/misc/utils.hpp"
#include "source/misc/universal_error.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"
#include "source/newtonian/test_1d/main_loop_1d.hpp"
#include "source/newtonian/one_dimensional/simple_cfl_1d.hpp"
#include "source/newtonian/one_dimensional/simple_extensive_updater_1d.hpp"
#include "source/newtonian/one_dimensional/simple_flux_calculator_1d.hpp"

// Riemann problem

using namespace std;
using namespace simulation1d;

class SimData
{
public:
  
  SimData(void):
    vertices_(linspace(0,1,100)),
    interpm_(),
    density_(1),
    pressure_(10,1,0.5),
    xvelocity_(0),
    yvelocity_(0),
    eos_(5./3.),
    rs_(eos_.getAdiabaticIndex()),
    vm_(true),
    bc_(),
    force_(),
    tsf_(0.3),
    fc_(rs_,interpm_,bc_),
    eu_(),
    cu_(),
    sim_(pg_,
	 vertices_,
	 density_,
	 pressure_,
	 xvelocity_,
	 yvelocity_,
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
  const vector<double> vertices_;
  PCM1D interpm_;
  const Uniform density_;
  const Step pressure_;
  const Uniform xvelocity_;
  const Uniform yvelocity_;
  const IdealGas eos_;
  const ERSIG rs_;
  const Lagrangian1D vm_;
  const RigidWall1D bc_;
  const ZeroForce1D force_;
  const SimpleCFL1D tsf_;
  const SimpleFluxCalculator1D fc_;
  const SimpleExtensiveUpdater1D eu_;
  const SimpleCellUpdater1D cu_;
  hdsim1D sim_;
};

namespace {
  void write_output(hdsim1D const& sim)
  {
    ofstream f("res.txt");
    const vector<Primitive>& cells = sim.getCells();
    f << cells.at(cells.size()/2).Pressure << endl;
    f << cells.at(cells.size()/2).Velocity.x << endl;
    f.close();
  }

#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
  __attribute__((noreturn))
#endif
void report_error(UniversalError const& eo)
{
  cout << "Caught universal error" << endl;
  cout << eo.GetErrorMessage() << endl;
  for(size_t i = 0;i<eo.GetFields().size();++i){
    cout << eo.GetFields()[i] << " = "
	 << eo.GetValues()[i] << endl;
  }
  throw;
}
}

int main(void)
{
  SimData sim_data;
  hdsim1D& sim = sim_data.getSim();

  try{
  main_loop(sim, 0.067, 1e6, 1,
	    "time.txt");
  }
  catch(UniversalError const& eo){
    report_error(eo);
  }
	
  write_output(sim);
	
  return 0;
}
