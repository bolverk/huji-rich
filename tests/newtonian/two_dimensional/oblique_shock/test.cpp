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
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/two_dimensional/condition_action_sequence.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/misc/vector_initialiser.hpp"

using namespace std;

namespace {

  vector<ComputationalCell> calc_init_cond
  (const Tessellation& tess)
  {
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i){
      res[i].density = 1;
      res[i].pressure = 1;
      res[i].velocity = Vector2D(4,0);
      const Vector2D r = tess.GetCellCM(static_cast<int>(i));
      res[i].stickers["wedge"] = r.y < 0.1*r.x-0.428;
    }
    return res;
  }

  vector<Vector2D> CustomPoints(int nx,int ny,double sidex,double sidey)
  {
    vector<Vector2D> res;
    double widthx = sidex/static_cast<double>(nx);
    double widthy = sidey/static_cast<double>(ny);
    Vector2D point;
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
	point.x = (static_cast<double>(i)+0.5)*widthx-sidex/2;
	point.y = (static_cast<double>(j)+0.5)*widthy-sidey/2;
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
    for(size_t i=0;i<inner_boundary_1.size();++i)
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
      pg_(),
      outer_(-width/2,width/2,width/2,-width/2),
      tess_(create_initial_points(width,np,inner_num),
	    outer_),
      eos_(adiabatic_index),
      rs_(),
      point_motion_(),
      force_(),
      tsf_(0.3),
      fc_(VectorInitialiser<pair<const ConditionActionSequence::Condition*,const ConditionActionSequence::Action*> >
	  (pair<const ConditionActionSequence::Condition*,const ConditionActionSequence::Action*>
	   (new IsBoundaryEdge, new FreeFlowFlux(rs_)))
	  (pair<const ConditionActionSequence::Condition*,const ConditionActionSequence::Action*>
	   (new RegularSpecialEdge("wedge"), new RigidWallFlux(rs_)))
	  (pair<const ConditionActionSequence::Condition*,const ConditionActionSequence::Action*>
	   (new IsBulkEdge, new RegularFlux(rs_)))()),
      eu_(),
      cu_(VectorInitialiser<pair<const SimpleCellUpdater::Condition*,const SimpleCellUpdater::Action*> >
	  (pair<const SimpleCellUpdater::Condition*,const SimpleCellUpdater::Action*>
	   (new HasSticker("wedge"),new SkipUpdate))()),
      sim_(tess_,
	   outer_,
	   pg_,
	   calc_init_cond(tess_),
	   eos_,
	   point_motion_,
	   force_,
	   tsf_,
	   fc_,
	   eu_,
	   cu_)
    {
      write_number(adiabatic_index,
		   "adiabatic_index.txt");
      write_number(4.0/sqrt(5./3.),
		   "mach_number.txt");
    }

    hdsim& getSim(void)
    {
      return sim_;
    }

  private:
    const SlabSymmetry pg_;
    SquareBox outer_;
    VoronoiMesh tess_;
    const IdealGas eos_;
    const Hllc rs_;
    Eulerian point_motion_;
    ZeroForce force_;
    const SimpleCFL tsf_;
    const ConditionActionSequence fc_;
    const SimpleExtensiveUpdater eu_;
    const SimpleCellUpdater cu_;
    hdsim sim_;
  };

  void main_loop(hdsim& sim)
  {
    SimpleCFL tsf(0.5);
    const int max_iter = 5e6;
    const double tf = 1;
    //    sim.SetEndTime(tf);
    while(tf>sim.getTime()){
      try{
	sim.TimeAdvance();
      }
      catch(UniversalError const& eo){
	DisplayError(eo);
      }

      if(sim.getCycle()>max_iter)
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


