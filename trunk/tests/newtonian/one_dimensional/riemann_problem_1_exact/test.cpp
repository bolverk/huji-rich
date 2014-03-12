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
    sim_(vertices_,
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
  hdsim1D sim_;
};

namespace {

void write_output(hdsim1D const& sim)
{
  ofstream f;
  f.open("res.txt");
  f << sim.GetCell(sim.GetCellNo()/2).Pressure << endl;
  f << sim.GetCell(sim.GetCellNo()/2).Velocity.x << endl;
  f.close();
}

void report_error(UniversalError const& eo)
{
  cout << "Caught universal error" << endl;
  cout << eo.GetErrorMessage() << endl;
  for(int i = 0;i<(int)eo.GetFields().size();++i){
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
