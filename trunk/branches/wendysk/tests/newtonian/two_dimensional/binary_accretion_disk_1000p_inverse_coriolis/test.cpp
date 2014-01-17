#define _USE_MATH_DEFINES
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
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
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
#include "source/newtonian/two_dimensional/pcm_scalar.hpp"
#include "source/misc/simple_io.hpp"

using namespace std;
using namespace simulation2d;

class AzimuthalVelocity: public SpatialDistribution
{
public:

  AzimuthalVelocity(Func1Var const& radial,
		    Vector2D const& center,
		    char comp):
    radial_(radial),
    center_(center),
    comp_(comp) {}

  double EvalAt(Vector2D const& p) const
  {
    const Vector2D rvec = p - center_;
    const double radius = abs(rvec);
    const double vq = radial_.eval(radius);
    if(comp_=='x')
      return -vq*rvec.y/radius;
    else if(comp_=='y')
      return vq*rvec.x/radius;
    else
      throw "Unknown component name "+comp_;
  }

private:
  Func1Var const& radial_;
  const Vector2D center_;
  const char comp_;
};

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

class RotatedKepler: public Func1Var
{
public:

  RotatedKepler(double donor_mass,
		double acceptor_mass,
		double separation):
    donor_mass_(donor_mass),
    acceptor_mass_(acceptor_mass),
    separation_(separation) {}

  double eval(double r) const
  {
    return r*(sqrt(acceptor_mass_/pow(r,3))-
	      sqrt((acceptor_mass_+donor_mass_)/pow(separation_,3)));
  }

private:
  const double donor_mass_;
  const double acceptor_mass_;
  const double separation_;
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
    pressure(1/adiabatic_index),
    donor_mass(1),
    donor_center(1,0),
    acceptor_mass(1),
    center_of_gravity((donor_center*acceptor_mass+
		       donor_mass*center)/
		      (donor_mass+acceptor_mass)),
    angular_velocity(sqrt(donor_mass+acceptor_mass)/
		     pow(abs(0.5*(donor_center-center)),3)),
    radial_velocity(donor_mass,
		    acceptor_mass,
		    abs(center-donor_center)) {}

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
  const RotatedKepler radial_velocity;
};

class AlignedRectangle
{
public:

  AlignedRectangle(Vector2D const& lower_left,
		   Vector2D const& upper_right):
    lower_left_(lower_left),
    upper_right_(upper_right) {}

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
	
	const double cell_size = sqrt(M_PI/mm.total_point_number)*mm.max_radius;
    const vector<Vector2D> v1 = CirclePointsRmax(mm.total_point_number,
						 mm.min_radius+cell_size,
						 mm.max_radius-cell_size,
						 mm.center.x,
						 mm.center.y);
    const int nr_outer = int(2*sqrt(mm.total_point_number*M_PI));
    const vector<Vector2D> v2 = Circle(nr_outer,
				       mm.max_radius-0.5*cell_size,
				       mm.center.x,
				       mm.center.y);
	 const vector<Vector2D> v3 = Circle(nr_outer,
				       mm.max_radius+0.5*cell_size,
				       mm.center.x,
				       mm.center.y);
    const int nr_inner = int(nr_outer*(mm.min_radius/mm.max_radius));
    const vector<Vector2D> v4 = Circle(nr_inner,
				       mm.min_radius+0.5*cell_size,
				       mm.center.x,
				       mm.center.y);
	const vector<Vector2D> v5 = Circle(nr_inner,
				       mm.min_radius-0.5*cell_size,
				       mm.center.x,
				       mm.center.y);
    return join(join(join(join(v1,v2),v3),v4),v5);    
  }

  class ConsecutiveSnapshots: public DiagnosticFunction
  {
  public:

    ConsecutiveSnapshots(double dt):
      next_time_(0),
      dt_(dt),
      counter_(0) {}

    void diagnose(hdsim const& sim)
    {
      write_number(sim.GetTime(),"time.txt");

      if(sim.GetTime()>next_time_){
	next_time_ += dt_;

	write_snapshot_to_hdf5(sim,"snapshot_"+int2str(counter_)+".h5");
	++counter_;
      }
    }

  private:
    double next_time_;
    const double dt_;
    int counter_;
  };

  void my_main_loop(hdsim& sim)
  {
    const double tf = 1;
    SafeTimeTermination term_cond(tf, 1e6);
    //    WriteTime diag("time.txt");
    ConsecutiveSnapshots diag(tf/100);
    try{
      main_loop(sim,
		term_cond,
		1,
		&diag);

      /*
		const double endtime=1;
		int counter=0;
		sim.SetEndTime(endtime);
		int dumpnumber=1;
		while(sim.GetTime()<endtime)
		{
			sim.TimeAdvance();
			if(counter%25==0)
			{
				cout<<"Time step number "<<counter<<" Time = "<<sim.GetTime()<<endl;
				write_snapshot_to_hdf5(sim,"c:\\sim_data\\jim"+int2str(dumpnumber)+".bin");
				++dumpnumber;
			}
			++counter;
		}
      */
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
    //point_motion_(pm_naive_),
	point_motion_(),
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
    zero_force_(),
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
	 AzimuthalVelocity(mm.radial_velocity,
			   mm.center,'x'),
	 AzimuthalVelocity(mm.radial_velocity,
			   mm.center,'y'),
	 eos_,
	 rs_,
	 &point_motion_,
	 //&zero_force_, 
	 &force_,
	 &outer_,
	 &hbc_,
	 false,
	 false) {}

  hdsim& getSim(void)
  {
    return sim_;
  }
   
private:
  VoronoiMesh tess_;
  PCM2D interp_;
  const IdealGas eos_;
  Lagrangian pm_naive_;
  //RoundCells point_motion_;
  Eulerian point_motion_;
  const Hllc rs_;
  CenterGravity acceptor_gravity_acceleration_;
  ConservativeForce acceptor_gravity_;
  CenterGravity donor_gravity_acceleration_;
  ConservativeForce donor_gravity_;
  Centripetal centripetal_acceleration_;
  ConservativeForce centripetal_force_;
  Coriolis coriolis_;
  SeveralSources force_;
  ZeroForce zero_force_;
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
    if (((dir_==in)&&(ScalarProd(ghost.Velocity,outward)>0))||
	((dir_==out)&&(ScalarProd(ghost.Velocity,outward)<0)))
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
    return FluxInBulk(n,p,left,right,face_velocity,rs);
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
  void report_error(UniversalError const& eo,
		    hdsim const& sim)
{
  cout << eo.GetErrorMessage() << endl;
  for(int i=0;i<(int)eo.GetFields().size();++i){
    cout << eo.GetFields()[i] << " = "
	 << eo.GetValues()[i] << endl;
    if(eo.GetFields()[i]=="cell index"){
      const int cell_index = (int)eo.GetValues()[i];
      cout << "position = " << sim.GetMeshPoint(cell_index).x << " , "
	   << sim.GetMeshPoint(cell_index).y << endl;
    }
  }
}
}

int main(void)
{
  SimData sim_data;
  hdsim& sim = sim_data.getSim();

 write_snapshot_to_hdf5(sim,"initial.h5");

  const MagicNumbers mm;
  Ratchet obc(Ratchet::in);
  for(int i=0;i<sim.GetCellNo();++i){
    const Vector2D pos = sim.GetMeshPoint(i);
    if(abs(pos)>mm.max_radius*1.00001)
      sim.CellsEvolve[i]=&obc;
    else if(abs(pos)<mm.min_radius*0.9999)
      sim.CellsEvolve[i]=&obc;
  }

  // Add Codl flows
  PCMScalar sinterp;
  sim.SetColdFlows(0.05,0.05);
  sim.setTracerInterpolation(&sinterp);

  try{
  my_main_loop(sim);
  }
  catch(UniversalError const& eo){
    report_error(eo,sim);
    throw;
  }
  
  write_snapshot_to_hdf5(sim,"final.h5");
}
