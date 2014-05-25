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
#include "source/newtonian/one_dimensional/periodic_1d.hpp"
#include "source/newtonian/one_dimensional/diagnostics_1d.hpp"
#include "source/misc/int2str.hpp"
#include "source/misc/utils.hpp"
#include "source/misc/universal_error.hpp"

using namespace std;

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
      throw "First argument in external file should either be e or l (short for eulerian and lagrangian vertex motion";

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
    density_(1,0.3,2,0.7,1),
    pressure_(1),
    xvelocity_(1),
    yvelocity_(0),
    lagrangian_(true),
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
	 eulerian_,
	 bc_) {}

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
  TwoSteps density_;
  Uniform pressure_;
  Uniform xvelocity_;
  Uniform yvelocity_;
  Hllc rs_;
  Eulerian1D eulerian_;
  Lagrangian1D lagrangian_;
  Periodic1D bc_;
  hdsim sim_;
};

vector<string> get_output_hydro_vars(void)
{
  vector<string> res;
  res.push_back("center");
  res.push_back("density");
  res.push_back("pressure");
  res.push_back("xvelocity");
  return res;
}

void write_snapshot(hdsim const& sim, int cycle)
{
  vector<string> output_vars = get_output_hydro_vars();
  for(int i=0;i<(int)output_vars.size();++i)
    write_cells_property(sim,output_vars[i],
			 output_vars[i]+"_list."+
			 int2str(cycle)+".txt",10); 
}

void main_loop(hdsim& sim, int tio)
{
  const double tf = 1;
  const int max_cycle = (int)1e6;
  Timer timer(0,tf/10);

  write_snapshot(sim,0);

  while(sim.GetTime()<tf){
    
    if(tio==0)
      sim.TimeAdvance();
    else
      sim.TimeAdvanceRK(tio);

    if(sim.GetCycle()>max_cycle){
      cout << "Max number of cycles exceeded" << endl;
      throw;
    }

    if(timer.isTime(sim.GetTime()))
      write_snapshot(sim,timer.getCycle());
  }
}

void write_output(hdsim const& sim)
{
  const int prec = 10;
  write_to_file(sim.GetTime(),"time.txt");
  write_cells_property(sim,"center","cell_centres.txt",prec);
  write_cells_property(sim,"density","densities.txt",prec);
  write_cells_property(sim,"pressure","pressures.txt",prec);
  write_cells_property(sim,"xvelocity","velocities.txt",prec);
}

int main(void)
{
  try{
  SimData sim_data(ExternalData("data_file.txt"));
  hdsim& sim = sim_data.getSim();
  //  sim.overrideCFL(0.01);
  
  main_loop(sim,sim_data.getTimeIntegrationOrder());

  write_output(sim);
  }
  catch(UniversalError const& eo){
    cout << eo.GetErrorMessage() << endl;
  }
  
  return 0;
}
