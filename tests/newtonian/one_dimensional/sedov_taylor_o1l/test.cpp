#include <iostream>
#include <cmath>
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/one_dimensional/lagrangian1d.hpp"
#include "source/newtonian/one_dimensional/rigid_wall_1d.hpp"
#include "source/newtonian/one_dimensional/cylindrical_complementary_1d.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"

namespace {

  vector<double> logspace
  (double base,
   double exp_min,
   double exp_max,
   size_t num)
  {
    vector<double> res(num,0);
    for(size_t i=0;i<num;++i)
      res.at(i) = pow
	(base,
	 exp_min+(exp_max-exp_min)/
	 static_cast<double>(num)*
	 static_cast<double>(i));
    return res;
  }

  class SimData
  {
  public:
    SimData(void):
      pg_(),
      interpm_(),
      eos_(5./3.),
      rs_(),
      vm_(false),
      bc_(),
      force_(),
      sim_
      (pg_,
       logspace(10,-2,0.5,200),
       interpm_,
       Uniform(1),
       Step(1e3,1e-9,1e-1),
       Uniform(0),
       Uniform(0),
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
    const SphericalSymmetry1D pg_;
    //    const SlabSymmetry1D pg_;
    const PCM1D interpm_;
    const IdealGas eos_;
    const Hllc rs_;
    const Lagrangian1D vm_;
    const RigidWall1D bc_;
    const CylindricalComplementary1D force_;
    hdsim1D sim_;
  };

  bool termination_condition
  (const hdsim1D& sim)
  {
    assert(sim.GetCellNo()>10);
    return sim.GetCell(sim.GetCellNo()-10).Density>1.2;
    //return sim.GetCycle()>1e5;
  }

  void my_main_loop(hdsim1D& sim)
  {
    while(!termination_condition(sim)){
      sim.TimeAdvance();
      //      std::cout << sim.GetTime() << " " << sim.GetCycle() << std::endl;
    }
    diagnostics1d::write_snapshot_to_hdf5(sim,"final.h5");
  }
}

int main(void) 
{
  SimData sim_data;
  hdsim1D& sim = sim_data.getSim();

  my_main_loop(sim);

  return 0;
}
