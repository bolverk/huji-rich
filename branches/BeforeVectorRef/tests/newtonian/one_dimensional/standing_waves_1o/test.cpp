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

using namespace std;
using namespace simulation1d;
using namespace diagnostics1d;

class PeriodicDriver: public ExternalForces1D
{
public:

  PeriodicDriver(double wavelength,
		 double amplitude,
		 double phase_velocity):
    k_(2*M_PI/wavelength),
    amp_(amplitude),
    v_(phase_velocity) {}

  Conserved calc
  (vector<double> const& vertices,
   vector<Primitive> const& cells,
   int point,
   double t,
   double /*dt*/) const
  {
    const double volume = vertices[point+1]-vertices[point];
    const double x = 0.5*(vertices[point+1]+vertices[point]);
    const double density = cells[point].Density;
    const double xvelocity =cells[point].Velocity.x;
    const double acceleration = amp_*sin(k_*x)*sin(k_*v_*t);
    const double xmom = density*acceleration;
    const double enr = density*acceleration*xvelocity;
    return -volume*Conserved(0,Vector2D(xmom,0),enr);
  }

private:
  const double k_;
  const double amp_;
  const double v_;
};

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
    sim_(linspace(0,width,30),
	 interpm_,
	 Uniform(read_number("mean_density.txt")),
	 Uniform(read_number("mean_pressure.txt")),
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
  PCM1D interpm_;
  const IdealGas eos_;
  const Hllc rs_;
  const Eulerian1D vm_;
  const RigidWall1D bc_;
  const PeriodicDriver force_;
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
