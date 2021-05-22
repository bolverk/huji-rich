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
#include "source/newtonian/three_dimensional/hdf5_diagnostics_3d.hpp"
#include <fstream>

using std::cout;
using std::endl;
using std::pair;
using std::ofstream;

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
  (double r_min,
   double r_max)
  {
    const double r = r_min*pow(r_max/r_min, random_float(0,1));
    const double q = acos(random_float(-1, 1));
    const double f = random_float(0, 2*M_PI);
    const double x = r*sin(q)*cos(f);
    const double y = r*sin(q)*sin(f);
    const double z = r*cos(q);
    return Vector3D(x,y,z);
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
   double r_min,
   double r_max)
  {
    vector<Vector3D> res(n);
    generate(res.begin(),
	     res.end(),
	     [&](){return generate_random_point(r_min, r_max);});
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
   const size_t n,
   const double inner_ratio=0.1)
  {
    const double outer_radius = abs(boundary.first - boundary.second);
    vector<Vector3D> raw = 
      generate_initial_grid(n,
			    outer_radius*inner_ratio,
			    outer_radius);
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
      const double r = abs(points.at(i));
      res.at(i).density = 1;
      res.at(i).pressure = r < 1e-2 ? 1e2 : 1e-9;
      res.at(i).velocity = Vector3D(0,0,0);
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
      boundary_(Vector3D(-1, -1, -1),
		Vector3D(1, 1, 1)),
      init_points_(create_grid(boundary_,
			      2000,
			      1e-3)),
      /*      tess_(boundary_.first,
	    boundary_.second),
      */
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
      tsn_({}, {}),
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
       tsn_) 
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
    TracerStickerNames tsn_;
    HDSim3D sim_;
  };

  double estimate_shock_radius
  (HDSim3D& sim,
   const double p_thres)
  {
    double res = 0;
    for(size_t i=0;i<sim.getCells().size();++i){
      if(sim.getCells().at(i).pressure>p_thres){
	const Vector3D r = sim.getTesselation().GetMeshPoint(i);
	const double radius = abs(r);
	res = std::max(res, radius);
      }
    }
    return res;
  }

  class TrackShockRadius
  {
  public:

    TrackShockRadius(const string& fname,
		     const double p_thres):
      fname_(fname), p_thres_(p_thres)  {}

    void operator()(HDSim3D& sim)
    {
      const double time = sim.getTime();
      const double radius =
	estimate_shock_radius(sim,
			      p_thres_);
      trajectory_.push_back({time, radius});
    }

    ~TrackShockRadius(void)
    {
      ofstream f(fname_.c_str());
      for(size_t i=0;i<trajectory_.size();++i)
	f << trajectory_.at(i).first << " "
	  << trajectory_.at(i).second << endl;
      f.close();
    }

  private:
    const string fname_;
    const double p_thres_;
    vector<pair<double, double> > trajectory_;
  };

  void main_loop(HDSim3D& sim)
  {
    const double tf = 5e-1;
    //    write_txt_snapshot(sim, "initial.txt");
    write_snapshot_to_hdf5(sim, "initial.h5");
    TrackShockRadius diag("shock_trajectory.txt",
			  1e-3);
    while(sim.getTime()<tf){
      sim.timeAdvance();
      diag(sim);
    }
    //    write_txt_snapshot(sim, "final.txt");
    write_snapshot_to_hdf5(sim, "final.h5");
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
