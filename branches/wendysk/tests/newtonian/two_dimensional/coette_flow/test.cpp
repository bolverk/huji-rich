#include <iostream>
#include <fstream>
#include <cmath>
#include "source/misc/utils.hpp"
#include "source/misc/func_1_var.hpp"
#include "source/tessellation/geometry.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/interpolations/pcm2d.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/two_dimensional/custom_evolutions/RigidBodyEvolve.hpp"
#include "source/newtonian/two_dimensional/custom_evolutions/ConstantPrimitiveEvolution.hpp"

using namespace std;
using namespace simulation2d;

class AlignedRectangle
{
public:

  AlignedRectangle(Vector2D const& lower_left,
		   Vector2D const& upper_right):
    lower_left_(lower_left),
    upper_right_(upper_right) {}

  bool contains(Vector2D const& point) const
  {
    return lower_left_.x < point.x and
      lower_left_.y < point.y and
      upper_right_.x > point.x and
      upper_right_.y > point.y;
  }

private:
  const Vector2D lower_left_;
  const Vector2D upper_right_;
};

namespace{

vector<Vector2D> radial_grid(vector<double> const& radii,
			     Vector2D const& center,
			     AlignedRectangle const& domain)
{
  vector<Vector2D> res;
  if(radii[0]==0)
    res.push_back(center);
  else{
    const double dq = M_PI/6.0;
    for(int i=0;i<6;++i)
      res.push_back(center+Vector2D(radii[0]*cos(i*dq),
				 radii[0]*sin(i*dq)));
  }    
  for(int i=1;i<(int)radii.size();++i){
    const double dr = radii[i]-radii[i-1];
    const double dq = dr/radii[i];
    const int nq = max(int(2*M_PI/dq),1);
    for(int j=0;j<nq;++j){
      const Vector2D candidate = 
	center + Vector2D(radii[i]*cos(2*M_PI*double(j)/double(nq)),
			  radii[i]*sin(2*M_PI*double(j)/double(nq)));
      if(domain.contains(candidate))
	res.push_back(candidate);
    }
  }
  return res;
}

 class RadialVelocity: public SpatialDistribution
 {
 public:

   RadialVelocity(Func1Var const& radial,
		  Vector2D const& center,
		  char comp):
     radial_(radial),
     center_(center),
     comp_(comp) {}

   double EvalAt(Vector2D const& p) const
   {
     const Vector2D rvec = p - center_;
     const double radius = abs(rvec);
     const double vr = radial_.eval(radius);
     if(comp_=='x')
       return vr*rvec.x/radius;
     else if(comp_=='y')
       return vr*rvec.y/radius;
     else
       throw "Unknown component name "+comp_;
   }

 private:
   Func1Var const& radial_;
   const Vector2D center_;
   const char comp_;
 };

 class BoundedPowerLaw: public Func1Var
 {
 public:

   BoundedPowerLaw(double prefactor, 
		   double power,
		   double upper_bound,
		   double lower_bound):
     prefactor_(prefactor),
     power_(power),
     upper_bound_(upper_bound),
     lower_bound_(lower_bound) {}

   double eval(double x) const
   {
     const double res = prefactor_*pow(x,power_);
     if(res>upper_bound_)
       return upper_bound_;
     else if(res<lower_bound_)
       return lower_bound_;
     else
       return res;
   }

 private:
   const double prefactor_;
   const double power_;
   const double upper_bound_;
   const double lower_bound_;
 };

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.1, 1e6);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      1,
	      &diag);
  }
}

class SimData
{
public:

  SimData(Vector2D const& lower_left = Vector2D(-0.5,-0.5),
	  Vector2D const& upper_right = Vector2D(0.5, 0.5),
	  vector<double> const& radii = linspace(0,0.45,50),
	  Vector2D const& center = Vector2D(0,0),
	  BoundedPowerLaw const& rad_vel_prof = BoundedPowerLaw(1e3,-0.5,1,1e-3)):
    tess_(),
    interp_(),
    eos_(5./3.),
    pm_naive_(),
    point_motion_(pm_naive_),
    rs_(),
    force_(),
    outer_(lower_left.x,
	   upper_right.y,
	   upper_right.x,
	   lower_left.y),
    hbc_(rs_),
    sim_(radial_grid(radii,
		     center,
		     AlignedRectangle(lower_left,
				      upper_right)),
	 &tess_,
	 &interp_,
	 Uniform2D(1),
	 Uniform2D(1),
	 RadialVelocity(rad_vel_prof,
			center,'x'),
	 RadialVelocity(rad_vel_prof,
			center,'y'),
	 eos_,
	 rs_,
	 &point_motion_,
	 &force_,
	 &outer_,
	 &hbc_) {};

  hdsim& getSim(void)
  {
    return sim_;
  }
	 
private:
  VoronoiMesh tess_;
  PCM2D interp_;
  const IdealGas eos_;
  Lagrangian pm_naive_;
  RoundCells point_motion_;
  const Hllc rs_;
  ZeroForce force_;
  const SquareBox outer_;
  const RigidWallHydro hbc_;
  hdsim sim_;
};

int main(void)
{
  SimData sim_data;
  hdsim& sim = sim_data.getSim();
  
  RigidBodyEvolve rigid;
  ConstantPrimitiveEvolution cpe;
  for(int i=0;i<sim.GetCellNo();++i){
    const Vector2D pos = sim.GetMeshPoint(i);
    if(abs(pos)>0.4)
      sim.CellsEvolve[i]=&cpe;
    else if(abs(pos)<0.1)
      sim.CellsEvolve[i]=&cpe;
  }

  my_main_loop(sim);
  
  write_snapshot_to_hdf5(sim,"final.h5");
}
