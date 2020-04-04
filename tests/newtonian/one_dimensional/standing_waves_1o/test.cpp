#include <iostream>
#include <cmath>
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/rigid_wall_1d.hpp"
#include "source/misc/simple_io.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/one_dimensional/source_term_1d.hpp"
#include "source/newtonian/test_1d/main_loop_1d.hpp"
#include "source/newtonian/one_dimensional/hdf5_diagnostics1d.hpp"
#include "source/newtonian/one_dimensional/simple_cfl_1d.hpp"
#include "source/newtonian/one_dimensional/simple_extensive_updater_1d.hpp"
#include "source/newtonian/one_dimensional/simple_flux_calculator_1d.hpp"

using namespace std;
using namespace simulation1d;
using namespace diagnostics1d;

namespace {
  class PeriodicDriver: public SourceTerm1D
  {
  public:

    PeriodicDriver(double wavelength,
		   double amplitude,
		   double phase_velocity):
      k_(2*M_PI/wavelength),
      amp_(amplitude),
      v_(phase_velocity) {}

    Extensive operator()
    (const SimulationState1D& state,
     size_t point,
     double t,
     double /*dt*/) const
    {
      const vector<double>& vertices = state.getVertices();
      const vector<ComputationalCell>& cells = state.getCells();
      const size_t index = static_cast<size_t>(point);
      const double volume = vertices[index+1]-vertices[index];
      const double x = 0.5*(vertices[index+1]+vertices[index]);
      const double density = cells[index].density;
      const double xvelocity =cells[index].velocity.x;
      const double acceleration = amp_*sin(k_*x)*sin(k_*v_*t);
      const double xmom = density*acceleration;
      const double enr = density*acceleration*xvelocity;

      Extensive res;
      res.mass = 0;
      res.momentum = -volume*xmom*Vector2D(1,0);
      res.energy = -volume*enr;
      res.tracers = vector<double>(cells.at(0).tracers.size(),0);
      return res;
    }

  private:
    const double k_;
    const double amp_;
    const double v_;
  };
}

class SimData
{
public:

  SimData(double width = 1):
    interpm_(),
    eos_(read_number("adiabatic_index.txt")),
    rs_(),
    vm_(),
    bc_(),
    force_(read_number("wavelength.txt"),
	   read_number("amplitude.txt"),
	   read_number("phase_velocity.txt")),
    tsf_(0.33333333333),
    fc_(rs_,interpm_,bc_),
    eu_(),
    cu_(),
    sim_(pg_,
	 SimulationState1D
	 (linspace(0,width,30),
	  Uniform(read_number("mean_density.txt")),
	  Uniform(read_number("mean_pressure.txt")),
	  Uniform(0),
	  Uniform(0),
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
  PCM1D interpm_;
  const IdealGas eos_;
  const Hllc rs_;
  const Eulerian1D vm_;
  const RigidWall1D bc_;
  const PeriodicDriver force_;
  const SimpleCFL1D tsf_;
  const SimpleFluxCalculator1D fc_;
  const SimpleExtensiveUpdater1D eu_;
  const SimpleCellUpdater1D cu_;
  hdsim1D sim_;
};

int main(void)
{
  try{
    SimData sim_data;
    hdsim1D& sim = sim_data.getSim();

    main_loop(sim, 1, 1e6, 1, 
	      "time.txt");

    write_snapshot_to_hdf5(sim, "final.h5");
  }
  catch(string eo){
    cout << eo << endl;
  }

  return 0;
}
