#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <boost/assign/list_of.hpp>
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
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/source_terms/ConservativeForce.hpp"
#include "source/newtonian/two_dimensional/source_terms/CenterGravity.hpp"
#include "source/newtonian/two_dimensional/source_terms/SeveralSources.hpp"

using namespace std;
using namespace simulation2d;

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

class PowerLaw: public Func1Var
{
public:
  PowerLaw(double prefactor,
	   double power):
    prefactor_(prefactor),
    power_(power) {}

  double eval(double x) const
  {
    return prefactor_*pow(x,power_);
  }

private:

  const double prefactor_;
  const double power_;
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

class MagicNumbers
{
public:

  MagicNumbers(void):
    total_point_number(1000),
    min_radius(0.05),
    max_radius(0.36),
    cell_size(sqrt(M_PI*(pow(max_radius,2)-pow(min_radius,2))/
		   total_point_number)),
    center(0,0),
    adiabatic_index(5./3.),
    density(1),
    pressure(1),
    donor_mass(1),
    donor_center(1,0),
    acceptor_mass(1),
    center_of_gravity((donor_center*acceptor_mass+
		       donor_mass*center)/
		      (donor_mass+acceptor_mass)),
    angular_velocity(sqrt(donor_mass+acceptor_mass)/
		     abs(0.5*(donor_center-center))),
    radial_velocity(sqrt(acceptor_mass),-0.5) {}

  const int total_point_number;
  const double min_radius;
  const double max_radius;
  const double cell_size;
  const Vector2D center;
  const double adiabatic_index;
  const double density;
  const double pressure;
  const double donor_mass;
  const Vector2D donor_center;
  const double acceptor_mass;
  const Vector2D center_of_gravity;
  const double angular_velocity;
  const PowerLaw radial_velocity;
};

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

  template<class T> vector<T> join(vector<T> const& v1,
				   vector<T> const& v2)
  {
    vector<T> res(v1.size()+v2.size());
    copy(v1.begin(),v1.end(),res.begin());
    copy(v2.begin(),v2.end(),res.begin()+v1.size());
    return res;
  }

  vector<Vector2D> grid_func(MagicNumbers const& mm=MagicNumbers())
  {
    const vector<Vector2D> v1 = CirclePointsRmax(mm.total_point_number,
						 mm.min_radius,
						 mm.max_radius,
						 mm.center.x,
						 mm.center.y);
    const double cell_size = sqrt(M_PI/mm.total_point_number)*mm.max_radius;
    const int nr_outer = int(2*sqrt(mm.total_point_number*M_PI));
    const vector<Vector2D> v2 = Circle(nr_outer,
				       mm.max_radius+0.5*cell_size,
				       mm.center.x,
				       mm.center.y);
    const int nr_inner = int(nr_outer*(mm.min_radius/mm.max_radius));
    const vector<Vector2D> v3 = Circle(nr_inner,
				       mm.min_radius-0.5*cell_size,
				       mm.center.x,
				       mm.center.y);
    return join(join(v1,v2),v3);    
  }

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(1, 1e6);
    WriteTime diag("time.txt");
    try{
      main_loop(sim,
		term_cond,
		1,
		&diag);
    }
    catch(UniversalError& eo){
      eo.AddEntry("time",sim.GetTime());
      throw;
    }
  }
}

class Centripetal: public Acceleration
{
public:

  Centripetal(double angular_velocity,
	      Vector2D const& center):
    omega_(angular_velocity),
    center_(center) {}

  Vector2D Calculate(Tessellation const* tess,
		     vector<Primitive> const& /*cells*/,
		     int point,
		     vector<Conserved> const& /*fluxes*/,
		     vector<Vector2D> const& /*point_velocity*/,
		     HydroBoundaryConditions const* /*hbc*/,
		     double /*time*/,
		     double /*dt*/)
  {
    const Vector2D rvec(tess->GetCellCM(point) - center_);
    return pow(omega_,2)*rvec;
  }

private:

  const double omega_;
  const Vector2D center_;
};

namespace {
  Vector2D zcross(Vector2D const& v)
  {
    return Vector2D(v.y,-v.x);
  }
}

class Coriolis: public SourceTerm
{
public:

  Coriolis(double angular_velocity):
    omega_(angular_velocity) {}

  Conserved Calculate(Tessellation const* tess,
		      vector<Primitive> const& cells,
		      int point,
		      vector<Conserved> const& /*fluxes*/,
		      vector<Vector2D> const& /*point_velocity*/,
		      HydroBoundaryConditions const* /*hbc*/,
		      vector<vector<double> > const& /*tracers*/,
		      vector<double>& /*dtracer*/,
		      double /*time*/,
		      double /*dt*/)
  {
    const Vector2D velocity = cells[point].Velocity;
    const Vector2D acceleration = 2*omega_*zcross(velocity);
    const double mass = tess->GetVolume(point)*cells[point].Density;
    const Vector2D force = mass*acceleration;
    return Conserved(0,force,0);
  }

private:
  const double omega_;
};

template<class T> vector<T> four_vector(T t1,
					T t2,
					T t3,
					T t4)
{
  vector<T> res;
  res.push_back(t1);
  res.push_back(t2);
  res.push_back(t3);
  res.push_back(t4);
  return res;
}

class SimData
{
public:

  SimData(MagicNumbers const& mm = MagicNumbers(),
	  Vector2D const& lower_left = Vector2D(-0.38,-0.38),
	  Vector2D const& upper_right = Vector2D(0.38, 0.38)):
    tess_(),
    interp_(),
    eos_(5./3.),
    pm_naive_(),
    point_motion_(pm_naive_),
    rs_(),
    acceptor_gravity_acceleration_
    (mm.acceptor_mass,
     0,mm.center),
    acceptor_gravity_(&acceptor_gravity_acceleration_),
    donor_gravity_acceleration_
    (mm.donor_mass,
     0,mm.donor_center),
    donor_gravity_(&donor_gravity_acceleration_),
    centripetal_acceleration_(mm.angular_velocity,
			      mm.center_of_gravity),
    centripetal_force_(&centripetal_acceleration_),
    coriolis_(mm.angular_velocity),
    force_(four_vector<SourceTerm*>(&donor_gravity_,
				    &acceptor_gravity_,
				    &centripetal_force_,
				    &coriolis_)),
    outer_(lower_left.x,
	   upper_right.y,
	   upper_right.x,
	   lower_left.y),
    hbc_(rs_),
    sim_(grid_func(),  
	 &tess_,
	 &interp_,
	 Uniform2D(mm.density),
	 Uniform2D(mm.pressure),
	 RadialVelocity(mm.radial_velocity,
			mm.center,'x'),
	 RadialVelocity(mm.radial_velocity,
			mm.center,'y'),
	 eos_,
	 rs_,
	 &point_motion_,
	 &force_,
	 &outer_,
	 &hbc_) {}

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
  CenterGravity acceptor_gravity_acceleration_;
  ConservativeForce acceptor_gravity_;
  CenterGravity donor_gravity_acceleration_;
  ConservativeForce donor_gravity_;
  Centripetal centripetal_acceleration_;
  ConservativeForce centripetal_force_;
  Coriolis coriolis_;
  SeveralSources force_;
  const SquareBox outer_;
  const RigidWallHydro hbc_;
  hdsim sim_;
};

class Ratchet: public CustomEvolution
{
public:

  enum DIRECTION {in, out};

  Ratchet(DIRECTION dir):
    dir_(dir) {}

  Conserved CalcFlux(Tessellation const* tess,
		     vector<Primitive> const& cells,
		     double /*dt*/,
		     SpatialReconstruction* /*interp*/,
		     Edge const& edge,
		     Vector2D const& face_velocity,
		     RiemannSolver const& rs,
		     int index,
		     HydroBoundaryConditions const* hbc,
		     double /*time*/,
		     vector<vector<double> > const& /*tracers*/)
  {
    int other;
    int my_index;
    if(edge.GetNeighbor(0)==index){
      my_index = 0;
      other = edge.GetNeighbor(1);
    }
    else{
      my_index = 1;
      other = edge.GetNeighbor(0);
    }
    if(hbc->IsGhostCell(other,tess))
      return Conserved();
    const Vector2D p = Parallel(edge);
    const Vector2D n = 
      tess->GetMeshPoint(edge.GetNeighbor(1))
      -tess->GetMeshPoint(edge.GetNeighbor(0));
    const Vector2D outward = 
      tess->GetMeshPoint(edge.GetNeighbor(1-my_index))
      -tess->GetMeshPoint(edge.GetNeighbor(my_index));
    Primitive ghost = cells[other];
    if (((dir_==in)and(ScalarProd(ghost.Velocity,outward)>0))or
	((dir_==out)and(ScalarProd(ghost.Velocity,outward)<0)))
      ghost.Velocity = Reflect(cells[other].Velocity, p);
    Primitive left, right;
    if(0==my_index){
      left = ghost;
      right = cells[other];
    }
    else{
      left = cells[other];
      right = ghost;
    }
    Conserved res = rs.Solve(left, right,
			     Projection(face_velocity,n));
    res.Momentum = res.Momentum.x*n/abs(n) + 
      res.Momentum.y*p/abs(p);
    return res;
  }

  Primitive UpdatePrimitive(vector<Conserved> const& /*intensives*/,
			    EquationOfState const* /*eos*/,
			    vector<Primitive>& cells,
			    int index)
  {
    return cells[index];
  }			    

private:
  const DIRECTION dir_;
};

namespace {
void report_error(UniversalError const& eo)
{
  cout << eo.GetErrorMessage() << endl;
  for(int i=0;i<(int)eo.GetFields().size();++i)
    cout << eo.GetFields()[i] << " = "
	 << eo.GetValues()[i] << endl;
}
}

int main(void)
{
  SimData sim_data;
  hdsim& sim = sim_data.getSim();

  const MagicNumbers mm;
  Ratchet obc(Ratchet::out);
  for(int i=0;i<sim.GetCellNo();++i){
    const Vector2D pos = sim.GetMeshPoint(i);
    if(abs(pos)>mm.max_radius-mm.cell_size)
      sim.CellsEvolve[i]=&obc;
    else if(abs(pos)<mm.min_radius+mm.cell_size)
      sim.CellsEvolve[i]=&obc;
  }

  try{
  my_main_loop(sim);
  }
  catch(UniversalError const& eo){
    report_error(eo);
    throw;
  }
  
  write_snapshot_to_hdf5(sim,"final.h5");
}
