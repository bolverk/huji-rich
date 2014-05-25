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
#include "source/newtonian/one_dimensional/diagnostics_1d.hpp"
#include "source/newtonian/one_dimensional/source_term_1d.hpp"

using namespace std;

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

  SimData(int res, double width = 1):
    interpm_(),
    eos_(read_number("adiabatic_index.txt")),
    rs_(),
    vm_(),
    bc_(),
    force_(read_number("wavelength.txt"),
	   read_number("amplitude.txt"),
	   read_number("phase_velocity.txt")),
    sim_(linspace(0,width,res),
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

namespace{
void main_loop(hdsim1D& sim)
{
  const int max_iter = (int)1e6;
  const double tf = 1;

  while(sim.GetTime()<tf){
    sim.TimeAdvance();

    write_number(sim.GetTime(),"time.txt");
    
    if(sim.GetCycle()>max_iter)
      throw "Error in main_loop: max number of iterations exceeded";
  }
}

void write_output(hdsim1D const& sim)
{
  const int prec = 14;
  write_cells_property(sim,"center","cell_centers.txt",prec);
  write_cells_property(sim,"density","densities.txt",prec);
  write_cells_property(sim,"pressure","pressures.txt",prec);
  write_cells_property(sim,"xvelocity","velocities.txt",prec);
}
}

int main(void)
{
  try{
    SimData sim_data(read_int("resolution.txt"));
  hdsim1D& sim = sim_data.getSim();

  main_loop(sim);

  write_output(sim);
  }
  catch(string eo){
    cout << eo << endl;
  }

  return 0;
}
