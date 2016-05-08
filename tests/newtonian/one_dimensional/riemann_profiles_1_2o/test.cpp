#include <iostream>
#include <fstream>
#include <cmath>
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/rigid_wall_1d.hpp"
#include "source/newtonian/one_dimensional/plm1d.hpp"
#include "source/newtonian/one_dimensional/eos_consistent1d.hpp"
#include "source/misc/universal_error.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/one_dimensional/zero_force_1d.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"

// Riemann problem

using namespace std;
using namespace interpolations1d;
using namespace diagnostics1d;

namespace {
double GetXVelocityAt(const hdsim1D& sim, double x)
{
  for(size_t i=1;i<static_cast<size_t>(sim.GetCellNo());i++){
    if(sim.GetCellCenter(i-1)<x&&sim.GetCellCenter(i)>x)
      return 0.5*(sim.GetCell(i-1).Velocity.x+
		  sim.GetCell(i).Velocity.x);
  }
  throw "x out of range in GetXVelocityAt";
}
}

class StopCond
{
public:

  StopCond(const hdsim1D& sim, double xl, double xr):
    sim_(sim), xl_(xl), xr_(xr), iter_(0), max_iter_(50000) {}

  bool TimeToStop(void)
  {
    const double tol = 0.01;

    if(iter_>max_iter_)
      throw "Maximum number of cycles exceede: Time step may approach zero";
    iter_++;

    const bool cond1 = abs(GetXVelocityAt(sim_,xl_))>tol;
    const bool cond2 = abs(GetXVelocityAt(sim_,xr_))>tol;
    const bool cond3 = (cond1||cond2);
    return cond3;
  }
    
private:

  const hdsim1D& sim_;
  const double xl_;
  const double xr_;
  int iter_;
  const int max_iter_;
};

namespace {
void ReportError(UniversalError& eo)
{
  cout << eo.GetErrorMessage() << endl;
  for(size_t i=0;i<eo.GetFields().size();i++){
    cout << eo.GetFields()[i] << " = "
	 << eo.GetValues()[i] << endl;
  }
}
}

class SimData
{
public:

  SimData(void):
    eos_(5./3.),
    plm_naive_(),
    interpm_(plm_naive_,eos_),
    rs_(),
    vm_(),
    bc_(),
    force_(),
    sim_(pg_,
	 linspace(0,1,read_int("resolution.txt")),
	 interpm_,
	 Uniform(1),
	 Step(10,1,0.5),
	 Uniform(0),
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
  const IdealGas eos_;
  PLM1D plm_naive_;
  EOSConsistent interpm_;
  const Hllc rs_;
  const Eulerian1D vm_;
  const RigidWall1D bc_;
  const ZeroForce1D force_;
  hdsim1D sim_;
};

namespace {
void main_loop(hdsim1D& sim)
{
  StopCond sc(sim,0.1,0.9);
  while(!sc.TimeToStop()){
    try{
      sim.TimeAdvance2();
    }
    catch(UniversalError& eo){
      eo.AddEntry("time",sim.GetTime());
      eo.AddEntry("cycle",sim.GetCycle());
      ReportError(eo);
      throw;
    }
  }
}
}

int main(void)
{
  SimData sim_data;
  hdsim1D& sim = sim_data.getSim();

  main_loop(sim);
 
  write_snapshot_to_hdf5(sim,"final.h5");  
	
  return 0;
}
