#include <iostream>
#include "source/newtonian/test_1d/acoustic.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/PeriodicBox.hpp"
#include "source/newtonian/two_dimensional/PeriodicHydro.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/linear_gauss.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/test_2d/profile_1d.hpp"
#include "source/newtonian/two_dimensional/lagrangian.hpp"
#include "source/newtonian/two_dimensional/round_cells.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/zero_force.hpp"
#include "source/newtonian/test_2d/square_grid.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/convergence/external_data.hpp"
#include "source/newtonian/two_dimensional/pcm2d.hpp"
#include "source/newtonian/two_dimensional/eos_consistent.hpp"
#include "source/newtonian/two_dimensional/eulerian.hpp"
#include "source/misc/simple_io.hpp"

using namespace std;

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

  SimData(ExternalData const& exd, int res):
    time_int_order_(exd.getTimeIntegrationOrder()),
    width_(1),
    init_points_(square_grid(width_,res)),
    outer_(0,width_,width_,0),
    eos_(5./3.),
    plm_naive_(outer_,&hbc_,true,false),
    plm_(plm_naive_,eos_),
    interpm_(choose_between<SpatialReconstruction>
	     (exd.getInterpolationFlag(),
	      plm_,pcm_)),
    init_cond_(1,3./5.,eos_,1e-6,width_),
    density_(init_cond_.getProfile("density")),
    pressure_(init_cond_.getProfile("pressure")),
    xvelocity_(init_cond_.getProfile("xvelocity")),
    yvelocity_(init_cond_.getProfile("yvelocity")),
    lagrangian_(lagrangian_naive_),
    hbc_(rs_),
    sim_(init_points_,
	 &tess_,
	 &interpm_,
	 density_,
	 pressure_,
	 xvelocity_,
	 yvelocity_,
	 &eos_,
	 &rs_,
	 &choose_between<PointMotion>
	 (exd.getVertexMotionFlag(),
	  lagrangian_,eulerian_),
	 &force_,
	 &outer_,
	 &hbc_) {}

  hdsim& getSim(void)
  {
    return sim_;
  }

  int getTimeIntegrationOrder(void)
  {
    return time_int_order_;
  }

private:

  const int time_int_order_;
  double width_;
  vector<Vector2D> init_points_;
  PeriodicBox outer_;
  VoronoiMesh tess_;
  IdealGas eos_;
  PCM2D pcm_;
  LinearGauss plm_naive_;
  EOSConsistent plm_;
  SpatialReconstruction& interpm_;
  AcousticInitCond init_cond_;
  Profile1D density_;
  Profile1D pressure_;
  Profile1D xvelocity_;
  Profile1D yvelocity_;
  Eulerian eulerian_;
  Lagrangian lagrangian_naive_;
  RoundCells lagrangian_;
  Hllc rs_;
  PeriodicHydro hbc_;
  ZeroForce force_;
  hdsim sim_;
};

void main_loop(hdsim& sim,
	       int tio)
{
  const int max_iter = (int)1e6;
  const double tf = 1;
  
  while(sim.GetTime()<tf){

    if(1==tio)
      sim.TimeAdvance();
    else if(2==tio)
      sim.TimeAdvance2Mid();
    else
      throw "Unknown time integration order index";

    if(sim.GetCycle()>max_iter){
      cout << "Max number of time advance cycles exceeded" << endl;
      throw;
    }
  }
}

int main(void)
{
  try{
  SimData sim_data
    (ExternalData("data_file.txt"),
     read_int("resolution.txt"));
  hdsim& sim = sim_data.getSim();

  write_x_plot(sim,"x_prof_initial.txt",14);

  main_loop(sim,
	    sim_data.getTimeIntegrationOrder());

  write_x_plot(sim,"x_prof_final.txt",14);
  write_number(sim.GetTime(),"time.txt");
  }
  catch(string const& eo){
    cout << eo << endl;
    throw;
  }

  return 0;
}
