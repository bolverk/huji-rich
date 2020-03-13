#include <iostream>
#include <fstream>
#include "source/newtonian/common/ersig.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/lagrangian1d.hpp"
#include "source/newtonian/one_dimensional/outflow1d.hpp"
#include "source/newtonian/one_dimensional/zero_force_1d.hpp"
#include "source/misc/int2str.hpp"
#include "source/misc/utils.hpp"
#include "source/misc/universal_error.hpp"
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
    pressure_(1),
    xvelocity_(-1,1,0.5),
    yvelocity_(0),
    eos_(5./3.),
    rs_(eos_.getAdiabaticIndex()),
    vm_(false),
    bc_(),
    force_(),
    cu_(),
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
     force_,
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
  const Uniform pressure_;
  const Step xvelocity_;
  const Uniform yvelocity_;
  const IdealGas eos_;
  const ERSIG rs_;
  const Lagrangian1D vm_;
  const Outflow bc_;
  const ZeroForce1D force_;
  const SimpleCellUpdater1D cu_;
  hdsim1D sim_;
};

namespace {

void write_output(hdsim1D const& sim, string const& fname)
{
  ofstream f(fname.c_str());
  f << sim.GetCell(static_cast<size_t>(sim.GetCellNo()/2)).Pressure << endl;
  f << sim.GetCell(static_cast<size_t>(sim.GetCellNo()/2)).Velocity.x << endl;
  f.close();
}
}

int main(void)
{
  SimData sim_data;
  hdsim1D& sim = sim_data.getSim();
  
  main_loop(sim, 0.22, 1e6, 1,
	    "time.txt");
	
  write_output(sim,"res.txt");

  return 0;
}
