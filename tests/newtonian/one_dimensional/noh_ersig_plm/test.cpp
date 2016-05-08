#include <iostream>
#include "source/newtonian/one_dimensional/rigid_wall_1d.hpp"
#include "source/newtonian/one_dimensional/outflow1d.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/one_dimensional/plm1d.hpp"
#include "source/newtonian/one_dimensional/eos_consistent1d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/common/ersig.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/one_dimensional/zero_force_1d.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"
#include "source/newtonian/test_1d/main_loop_1d.hpp"

using namespace std;
using namespace interpolations1d;
using namespace simulation1d;
using namespace diagnostics1d;

namespace {
  class DifferentBC: public BoundaryConditions1D
  {
  public:

    DifferentBC(BoundaryConditions1D const& left,
		BoundaryConditions1D const& right):
      left_(left), right_(right) {}

    Conserved CalcFlux
    (vector<double> const& edges,
     vector<Primitive> const& cells,
     RiemannSolver const& rs,
     vector<double> const& vertex_velocity,
     int i) const
    {
      if(0==i)
	return left_.CalcFlux(edges,cells,rs,vertex_velocity,i);
      else if(static_cast<int>(edges.size())-1==i)
	return right_.CalcFlux(edges,cells,rs,vertex_velocity,i);
      else
	throw "Inside bulk of grid";
    }

  private:

    BoundaryConditions1D const& left_;
    BoundaryConditions1D const& right_;
  };
}

class SimData
{
public:

  SimData(void):
    pg_(),
    edges_(linspace(0,1,100)),
    pcm_(),
    plm_naive_(),
    density_(1),
    pressure_(1e-6),
    xvelocity_(-1),
    yvelocity_(0),
    eos_(5./3.),
    interpm_(plm_naive_,eos_),
    rs_(eos_.getAdiabaticIndex()),
    vm_(),
    left_bc_(),
    right_bc_(),
    bc_(left_bc_,right_bc_),
    force_(),
    sim_
    (pg_,
     edges_,
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
  const vector<double> edges_;
  const PCM1D pcm_;
  const PLM1D plm_naive_;
  const Uniform density_;
  const Uniform pressure_;
  const Uniform xvelocity_;
  const Uniform yvelocity_;
  const IdealGas eos_;
  const EOSConsistent interpm_;
  const ERSIG rs_;
  const Eulerian1D vm_;
  const RigidWall1D left_bc_;
  const Outflow right_bc_;
  const DifferentBC bc_;
  const ZeroForce1D force_;
  hdsim1D sim_;
};

int main(void)
{
  SimData sim_data;
  hdsim1D& sim = sim_data.getSim();

  main_loop(sim, 2, 1e6, 2, "time.txt");

  write_snapshot_to_hdf5(sim,"final.h5");

  return 0;
}
