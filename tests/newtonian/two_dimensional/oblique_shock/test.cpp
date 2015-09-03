#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "source/misc/int2str.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/tessellation/geometry.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/custom_evolutions/RigidBodyEvolve.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/FreeFlow.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/CustomOuter.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/InFlow.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"

using namespace std;

namespace {
  vector<Vector2D> CustomPoints(int nx,int ny,double sidex,double sidey)
  {
    vector<Vector2D> res;
    double widthx = sidex/(double)nx;
    double widthy = sidey/(double)ny;
    Vector2D point;
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
	point.x = ((double)i+0.5)*widthx-sidex/2;
	point.y = ((double)j+0.5)*widthy-sidey/2;
	if(point.y>0.1*point.x-0.428)
	  res.push_back(point);
      }
    }
    return res;
  }

  vector<Vector2D> Line(int PointNum,double xmin,double xmax,double ymin,double ymax)
  {
    double dy=ymax-ymin,dx=xmax-xmin;
    double angle=atan2(dy,dx);
    double length=sqrt(dy*dy+dx*dx);
    double dl=length/PointNum;

    vector<Vector2D> res;
    Vector2D temp;

    for(int i=0;i<PointNum;++i)
      {
	temp.Set(xmin+dl*(i+0.5)*cos(angle),ymin+dl*(i+0.5)*sin(angle));
	res.push_back(temp);
      }
    return res;
  }

  class HBCData
  {
  public:

    HBCData(Hllc const& rs,
	    Primitive const& hv):
      f_hbc_(rs),
      r_hbc_(rs),
      i_hbc_(hv,rs),
      hbc_(i_hbc_,f_hbc_,r_hbc_,f_hbc_)
    {
      write_number(hv.Velocity.x/hv.SoundSpeed,
		   "mach_number.txt");
    }

    CustomOuter& getHBC(void)
    {
      return hbc_;
    }

  private:
    FreeFlow f_hbc_;
    RigidWallHydro r_hbc_;
    InFlow i_hbc_;
    CustomOuter hbc_;
  };

  vector<Vector2D> create_initial_points
  (double width, int np, int inner_num)
  {
    vector<Vector2D> init_points = CustomPoints(np,np,width,width);
    vector<Vector2D> inner_boundary_1 = 
      Line(inner_num,
	   -width/2+0.01,
	   width/2,
	   -width/2,
	   width/2-0.9);
    const double alpha = atan((inner_boundary_1[0].y-inner_boundary_1[1].y)/
			      (inner_boundary_1[0].x-inner_boundary_1[1].x));
    write_number(alpha,"wedge_angle.txt");
    vector<Vector2D> inner_boundary_2(inner_boundary_1.size());
    const double a = 0.01;
    for(int i=0;i<(int)inner_boundary_1.size();++i)
      inner_boundary_2[i] = Vector2D(inner_boundary_1[i].x-a*sin(alpha),
				     inner_boundary_1[i].y+a*cos(alpha));
    init_points.insert
      (init_points.begin(), inner_boundary_2.begin(), inner_boundary_2.end());
    init_points.insert
      (init_points.begin(), inner_boundary_1.begin(), inner_boundary_1.end());
    return init_points;
  }

  class SimData
  {
  public:

    SimData(double adiabatic_index = 5./3.,
	    double width=1,
	    int np=30,
	    int inner_num=60):
      eos_(adiabatic_index),
      rs_(),
      point_motion_(),
      tess_(),
      hbc_data_(rs_,CalcPrimitive(1,1,Vector2D(4,0),eos_)),
      outer_(-width/2,width/2,width/2,-width/2),
      interp_(eos_,outer_,hbc_data_.getHBC(),
	      true, false),
      force_(),
      rigid_(),
      sim_(create_initial_points(width,np,inner_num),
	   tess_,
	   interp_,
	   Uniform2D(1),
	   Uniform2D(1),
	   Uniform2D(4),
	   Uniform2D(0),
	   eos_,
	   rs_,
	   point_motion_,
	   force_,
	   outer_,
	   hbc_data_.getHBC())
    {
      sim_.custom_evolution_manager.addCustomEvolution(&rigid_);
      for(size_t i=0;i<(size_t)inner_num;++i)
	sim_.custom_evolution_indices[i] = 1;

      write_number(adiabatic_index,
		   "adiabatic_index.txt");
    }

    hdsim& getSim(void)
    {
      return sim_;
    }

  private:
    const IdealGas eos_;
    const Hllc rs_;
    Eulerian point_motion_;
    VoronoiMesh tess_;
    HBCData hbc_data_;
    SquareBox outer_;
    LinearGaussConsistent interp_;
    ZeroForce force_;
    RigidBodyEvolve rigid_;
    hdsim sim_;
  };

  void main_loop(hdsim& sim)
  {
    SimpleCFL tsf(0.5);
    sim.setTimeStepFunction(tsf);
    const int max_iter = 5e6;
    const double tf = 1;
    //    sim.SetEndTime(tf);
    while(tf>sim.GetTime()){
      try{
	sim.TimeAdvance2Mid();
      }
      catch(UniversalError const& eo){
	DisplayError(eo);
      }

      if(sim.GetCycle()>max_iter)
	throw UniversalError
	  ("Maximum number of iterations exceeded in main loop");
    }
  }
}

int main(void)
{
  SimData sim_data;

  hdsim& sim = sim_data.getSim();

  main_loop(sim);

  write_snapshot_to_hdf5(sim, "final.h5");

  return 0;
}


