#include <iostream>
#include "source/newtonian/one_dimensional/rigid_wall_1d.hpp"
#include "source/newtonian/one_dimensional/outflow1d.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/common/ersig.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/one_dimensional/zero_force_1d.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"
#include "source/newtonian/test_1d/main_loop_1d.hpp"
#include "source/newtonian/one_dimensional/simple_cfl_1d.hpp"
#include "source/newtonian/one_dimensional/simple_extensive_updater_1d.hpp"
#include "source/newtonian/one_dimensional/simple_flux_calculator_1d.hpp"

using namespace std;
using namespace simulation1d;
using namespace diagnostics1d;

namespace {
  class DifferentBC: public BoundaryConditions1D
  {
  public:

    DifferentBC(BoundaryConditions1D const& left,
		BoundaryConditions1D const& right):
      left_(left), right_(right) {}

    Extensive operator()
    (const SimulationState1D& ss,
     const EquationOfState& eos,
     const RiemannSolver& rs,
     const vector<double>& vertex_velocity,
     const size_t i) const
    {
      const vector<double>& vertices = ss.getVertices();
      if(0==i)
	return left_(ss,eos,rs,vertex_velocity,i);
      else if(vertices.size()-1==i)
	return right_(ss,eos,rs,vertex_velocity,i);
      else
	throw "Inside bulk of grid";
    }
    
  private:
    
    BoundaryConditions1D const& left_;
    BoundaryConditions1D const& right_;
  };
  
  class SimData
  {
  public:
    
    SimData(void):
      edges_(linspace(0,1,1000)),
      interpm_(),
      density_(1),
      pressure_(1e-6),
      xvelocity_(-1),
      yvelocity_(0),
      eos_(5./3.),
      rs_(eos_.getAdiabaticIndex()),
      vm_(),
      left_bc_(),
      right_bc_(),
      bc_(left_bc_,right_bc_),
      force_(),
      tsf_(0.3),
      fc_(rs_,interpm_,bc_),
      eu_(),
      cu_(),
      sim_(pg_,
	   SimulationState1D
	   (edges_,
	    density_,
	    pressure_,
	    xvelocity_,
	    yvelocity_,
	    vector<pair<string, const SpatialDistribution1D*> >(),
	    vector<pair<string, const BoolSpatialDistribution*> >()),
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
    const vector<double> edges_;
    const PCM1D interpm_;
    const Uniform density_;
    const Uniform pressure_;
    const Uniform xvelocity_;
    const Uniform yvelocity_;
    const IdealGas eos_;
    const ERSIG rs_;
    const Eulerian1D vm_;
    const RigidWall1D left_bc_;
    const Outflow right_bc_;
    const DifferentBC bc_;
    const ZeroForce1D force_;
    const SimpleCFL1D tsf_;
    const SimpleFluxCalculator1D fc_;
    const SimpleExtensiveUpdater1D eu_;
    const SimpleCellUpdater1D cu_;
    hdsim1D sim_;
  };
}

int main(void)
{
  SimData sim_data;
  hdsim1D& sim = sim_data.getSim();

  main_loop(sim, 2, 1e6, 2, "time.txt");

  write_snapshot_to_hdf5(sim,"final.h5");

  return 0;
}
