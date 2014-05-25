#include <iostream>
#include <fstream>
#include <cmath>
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/one_dimensional/hdsim.hpp"
#include "source/newtonian/one_dimensional/pcm1d.hpp"
#include "source/newtonian/one_dimensional/plm1d.hpp"
#include "source/newtonian/one_dimensional/eos_consistent.hpp"
#include "source/newtonian/one_dimensional/spatial_distribution1d.hpp"
#include "source/newtonian/one_dimensional/eulerian1d.hpp"
#include "source/newtonian/one_dimensional/lagrangian1d.hpp"
#include "source/newtonian/one_dimensional/rigid_wall_1d.hpp"
#include "source/newtonian/one_dimensional/outflow1d.hpp"
#include "source/newtonian/one_dimensional/diagnostics_1d.hpp"
#include "source/misc/int2str.hpp"
#include "source/misc/utils.hpp"

// Riemann problem

using namespace std;

void WriteTime(hdsim const& sim, string fname)
{
  ofstream f;
  f.open(fname.c_str());
  f << sim.GetTime() << endl;
  f.close();
}

double GetXVelocityAt(const hdsim* sim, double x)
{
  for(int i=1;i<sim->GetCellNo();i++){
    if(sim->GetCellCenter(i-1)<x&&sim->GetCellCenter(i)>x)
      return 0.5*(sim->GetCell(i-1).Velocity.x+
		  sim->GetCell(i).Velocity.x);
  }
  throw "x out of range in GetXVelocityAt";
}

class StopCond
{
public:

  StopCond(const hdsim* sim, double xl, double xr):
    sim_(sim), xl_(xl), xr_(xr), iter_(0), max_iter_(50000) {}

  bool TimeToStop(void)
  {
    const double tol = 0.01;

    if(iter_>max_iter_)
      throw "Maximum number of cycles exceede: Time step may approach zero";
    iter_++;

    bool cond1 = abs(GetXVelocityAt(sim_,xl_))>tol;
    bool cond2 = abs(GetXVelocityAt(sim_,xr_))>tol;
    bool cond3 = (cond1||cond2);
    return cond3;
  }
    
private:

  const hdsim* sim_;
  double xl_;
  double xr_;
  int iter_;
  int max_iter_;
};

class ExternalData
{
public:

  ExternalData(string const& fname)
  {
    char buf_c;
    int buf_i;
    ifstream f(fname.c_str());

    f >> buf_c;
    if('e'==buf_c)
      vertex_motion_flag_ = false;
    else if('l'==buf_c)
      vertex_motion_flag_ = true;
    else
      throw "First argument in external file should either be e or l (short for eulerian and lagrangian vertex motion)";

    f >> buf_c;
    if('f'==buf_c)
      interpolation_flag_ = false;
    else if('s'==buf_c)
      interpolation_flag_ = true;
    else
      throw "Second argument in external file should either be f or s (short for first or second order";

    f >> buf_i;
    time_int_ord_ = buf_i;

    f >> buf_i;
    if(buf_i>=10)
      point_number_ = buf_i;
    else
      throw "Number of points must be larger than 10";

    f.close();
  }

  bool getVertexMotionFlag(void) const
  {
    return vertex_motion_flag_;
  }

  bool getInterpolationFlag(void) const
  {
    return interpolation_flag_;
  }

  int getTimeIntegrationOrder(void) const
  {
    return time_int_ord_;
  }

  int getPointNumber(void) const
  {
    return point_number_;
  }

private:
  bool vertex_motion_flag_;
  bool interpolation_flag_;
  int time_int_ord_;
  int point_number_;
};

template<class T> T& choose_between(bool flag,
				    T& tres,
				    T& fres)
{
  if(flag)
    return tres;
  else
    return fres;
}

class SimData
{
public:

  SimData(ExternalData const& exd):
    time_int_ord_(exd.getTimeIntegrationOrder()),
    vertices_(linspace(0,1,exd.getPointNumber())),
    eos_(5./3.),
    plm_(plm_naive_,eos_),
    density_(1),
    pressure_(1),
    xvelocity_(-1,1,0.5),
    yvelocity_(0),
    lagrangian_(false),
    sim_(vertices_,
	 choose_between<SpatialReconstruction>
	 (exd.getInterpolationFlag(),
	  plm_,pcm_),
	 density_,
	 pressure_,
	 xvelocity_,
	 yvelocity_,
	 eos_,
	 rs_,
	 choose_between<VertexMotion>
	 (exd.getVertexMotionFlag(),
	  lagrangian_,eulerian_),
	 bc_)
  {}

  int getTimeIntegrationOrder(void)
  {
    return time_int_ord_;
  }

  hdsim& getSim(void)
  {
    return sim_;
  }

private:
  int time_int_ord_;
  vector<double> vertices_;
  IdealGas eos_;
  PCM1D pcm_;
  PLM1D plm_naive_;
  EOSConsistent plm_;
  Uniform density_;
  Uniform pressure_;
  Step xvelocity_;
  Uniform yvelocity_;
  Hllc rs_;
  Eulerian1D eulerian_;
  Lagrangian1D lagrangian_;
  Outflow bc_;
  hdsim sim_;
};

void main_loop(hdsim& sim, int tio)
{
  const double tf = 0.05;
  const int max_cycle = (int)1e6;
  while(tf>sim.GetTime()){

    if(0==tio)
      sim.TimeAdvance();
    else{
      sim.TimeAdvanceRK(tio);
    }

    if(sim.GetCycle()>max_cycle){
      cout << "Max number of cycles exceeded" << endl;
      throw;
    }
  }
}

void write_output(hdsim const& sim)
{
  write_to_file(sim.GetTime(),"time.txt");
  write_cells_property(sim,"center","cell_centres.txt");
  write_cells_property(sim,"density","densities.txt");
  write_cells_property(sim,"pressure","pressures.txt");
  write_cells_property(sim,"xvelocity","velocities.txt");
}

int main(void)
{
  SimData sim_data(ExternalData("data_file.txt"));
  hdsim& sim = sim_data.getSim();
			
  main_loop(sim,sim_data.getTimeIntegrationOrder());

  write_output(sim);
  
  return 0;
}
