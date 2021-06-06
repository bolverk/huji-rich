#include <iostream>
#include <algorithm>
#include "source/3D/GeometryCommon/Vector3D.hpp"
#include "source/3D/GeometryCommon/Voronoi3D.hpp"
#include "source/newtonian/three_dimensional/computational_cell.hpp"
#include "source/newtonian/three_dimensional/Lagrangian3D.hpp"
#include "source/newtonian/three_dimensional/RoundCells3D.hpp"
#include "source/newtonian/three_dimensional/SourceTerm3D.hpp"
#include "source/newtonian/three_dimensional/CourantFriedrichsLewy.hpp"
#include "source/newtonian/three_dimensional/ConditionActionFlux1.hpp"
#include "source/newtonian/three_dimensional/Hllc3D.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/three_dimensional/PCM3D.hpp"
#include "source/newtonian/three_dimensional/default_cell_updater.hpp"
#include "source/newtonian/three_dimensional/default_extensive_updater.hpp"
#include "source/newtonian/three_dimensional/hdsim_3d.hpp"
//#include "source/newtonian/three_dimensional/hdf5_diagnostics_3d.hpp"
#include "source/3D/GeometryCommon/hdf_write.hpp"
#include <fstream>

using std::cout;
using std::endl;
using std::pair;
using std::for_each;

namespace {

  double random_float
  (double low,
   double high)
  {
    const double x = static_cast<double>(rand())/
      static_cast<double>(RAND_MAX);
    return low+x*(high-low);
  }

  Vector3D generate_random_point
  (const pair<Vector3D, Vector3D>& opposite_corners)
  {
    auto comp = {&Vector3D::x,
		 &Vector3D::y,
		 &Vector3D::z};
    Vector3D res;
    for_each(comp.begin(),
	     comp.end(),
	     [&res, &opposite_corners](auto c)
	     {
	       res.*c = random_float
		 (opposite_corners.first.*c,
		  opposite_corners.second.*c);
	     });
    return res;
  }

  bool is_inside(const Vector3D& point,
		 const pair<Vector3D, Vector3D>& boundary)
  {
    return point.x > boundary.first.x &&
      point.y > boundary.first.y &&
      point.z > boundary.first.z &&
      point.x < boundary.second.x &&
      point.y < boundary.second.y &&
      point.z < boundary.second.z;
  }

  vector<Vector3D> generate_initial_grid
  (size_t n,
   const pair<Vector3D, Vector3D>& opposite_corners)
  {
    vector<Vector3D> res(n);
    generate(res.begin(),
	     res.end(),
	     [&](){return generate_random_point(opposite_corners);});
    return res;
  }

  vector<Vector3D> lloyd_iteration
  (const vector<Vector3D>& points,
   const pair<Vector3D, Vector3D>& boundary,
   const double weight=0.1)
  {
    Voronoi3D tess(boundary.first,
		   boundary.second);
    tess.Build(points);
    vector<Vector3D> res(points.size());
    for(size_t i=0;i<points.size();++i)
      res.at(i) = tess.GetMeshPoint(i) + 
	weight*(tess.GetCellCM(i)- tess.GetMeshPoint(i));
    return res;
  }

  vector<Vector3D> create_grid
  (const pair<Vector3D, Vector3D>& boundary,
   const size_t n)
  {
    vector<Vector3D> raw = 
      generate_initial_grid
      (n,
       boundary);
    vector<Vector3D> res;
    for(size_t i=0;i<raw.size();++i){
      if(is_inside(raw.at(i), boundary))
	res.push_back(raw.at(i));
    }
    for(size_t i=0;i<50;++i)
      res = lloyd_iteration(res, boundary);
    return res;
  }

  vector<ComputationalCell3D> calc_init_cond
  (const vector<Vector3D>& points,
   const EquationOfState& eos)
  {
    vector<ComputationalCell3D> res(points.size());
    for(size_t i=0;i<res.size();++i){
      const Vector3D& r = points.at(i);
      const bool flag = abs(r)<0.3;
      res.at(i).density = 1;
      res.at(i).pressure = flag ? 2 : 1;
      res.at(i).velocity = flag ? Vector3D(0.1,0.2,-0.3) : Vector3D(0,0,0);
      res.at(i).internal_energy =
	eos.dp2e(res.at(i).density,
		 res.at(i).pressure);
    }
    return res;
  }

  class FluxCalculatorData
  {
  public:

    FluxCalculatorData(void):
      rs_(),
      cond_bulk_(),
      action_bulk_(rs_),
      cond_boundary_(),
      action_boundary_(rs_),
      ghost_(),
      interp_(ghost_),
      fc_({{&cond_bulk_, &action_bulk_},
	    {&cond_boundary_, &action_boundary_}},
	interp_) {}

    ConditionActionFlux1& getFluxCalculator(void)
    {
      return fc_;
    }

  private:
    Hllc3D rs_;
    IsBulkFace3D cond_bulk_;
    RegularFlux3D action_bulk_;
    IsBoundaryFace3D cond_boundary_;
    RigidWallFlux3D action_boundary_;
    RigidWallGenerator3D ghost_;
    PCM3D interp_;
    ConditionActionFlux1 fc_;
  };

  class TessellationData
  {
  public:

    TessellationData
    (const pair<Vector3D, Vector3D>& boundary,
     const vector<Vector3D>& points):
      tess(boundary.first,
	   boundary.second)
    {
      tess.Build(points);
    }

    Voronoi3D tess;
  };

class SimData
  {
  public:

    SimData(void):
      boundary_
      (Vector3D(-1, -1, -1),
       Vector3D(1, 1, 1)),
      init_points_
      (create_grid
       (boundary_,
	2000)),
      td_(boundary_, init_points_),
      eos_(5./3.),
      bpm_(),
      point_motion_(bpm_,
		    eos_,
		    boundary_.first,
		    boundary_.second),
      source_term_(),
      tsf_(0.3, 0.3, source_term_),
      fcd_(),
      cu_(),
      eu_(),
      sim_
      (td_.tess,
       calc_init_cond(init_points_, eos_),
       eos_,
       point_motion_,
       tsf_,
       fcd_.getFluxCalculator(),
       cu_,
       eu_,
       source_term_,
      {vector<string>(), vector<string>()}) 
    {}

    HDSim3D& getSim(void)
    {
      return sim_;
    }

  private:
    const pair<Vector3D, Vector3D> boundary_;
    const vector<Vector3D> init_points_;
    //    Voronoi3D tess_;
    TessellationData td_;
    const IdealGas eos_;
    Lagrangian3D bpm_;
    RoundCells3D point_motion_;
    ZeroForce3D source_term_;
    CourantFriedrichsLewy tsf_;
    FluxCalculatorData fcd_;
    DefaultCellUpdater cu_;
    DefaultExtensiveUpdater eu_;
    HDSim3D sim_;
  };

  void main_loop(HDSim3D& sim)
  {
    const double tf = 0.05;
    // write_snapshot_to_hdf5(sim, "initial.h5");
    WriteSnapshot3D(sim, "initial.h5");
    while(sim.getTime()<tf){
      sim.timeAdvance();
    }
    //   write_snapshot_to_hdf5(sim, "final.h5");
    WriteSnapshot3D(sim, "final.h5");
  }
}

int main(void)
{
  try{
    main_loop(SimData().getSim());
  }
  catch(const UniversalError& eo){
    reportError(eo);
  }
  return 0;
}
