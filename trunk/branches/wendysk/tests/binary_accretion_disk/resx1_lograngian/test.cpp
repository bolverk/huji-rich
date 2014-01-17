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
#include "source/newtonian/two_dimensional/interpolations/linear_gauss.hpp"
#include "source/newtonian/two_dimensional/interpolations/eos_consistent.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
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
#include "source/newtonian/two_dimensional/RemovalStrategy.hpp"
#include "source/newtonian/two_dimensional/linear_gauss_scalar.hpp"
#include "source/newtonian/two_dimensional/RefineStrategy.hpp"
#include "source/newtonian/two_dimensional/wendysk/polar2cart.hpp"
#include "source/newtonian/two_dimensional/wendysk/disk_amr.hpp"
#include "source/newtonian/two_dimensional/wendysk/centripetal.hpp"
#include "source/newtonian/two_dimensional/wendysk/coriolis.hpp"
#include "source/newtonian/two_dimensional/wendysk/ratchet.hpp"
#include "source/newtonian/two_dimensional/wendysk/const_primitive.hpp"
#include "source/newtonian/two_dimensional/wendysk/edge_repellant.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "source/newtonian/two_dimensional/Reset.hpp"

using namespace std;
using namespace simulation2d;

class MagicNumbers
{
public:

	MagicNumbers(void):
	  center(0,0),
		  adiabatic_index(1.2),
		  density(1e-5),
		  pressure(1e-4/adiabatic_index),
		  donor_mass(0.5),
		  donor_center(1,0),
		  acceptor_mass(0.5),
		  center_of_gravity((donor_center*acceptor_mass+
		  donor_mass*center)/
		  (donor_mass+acceptor_mass)),
		  angular_velocity(sqrt(donor_mass+acceptor_mass)/
		  pow(abs((donor_center-center)),3)),
		  orifice(0.005),
		  wind_pressure(1e-2/adiabatic_index),
		  wind_velocity(-1e-2,0),
		  wind_density(1) {}

	  const Vector2D center;
	  const double adiabatic_index;
	  const double density;
	  const double pressure;
	  const double donor_mass;
	  const Vector2D donor_center;
	  const double acceptor_mass;
	  const Vector2D center_of_gravity;
	  const double angular_velocity;
	  const double orifice;
	  const double wind_pressure;
	  const Vector2D wind_velocity;
	  const double wind_density;
};

class GridData
{
public:

	GridData(void):
	  lower_left(-0.5,-0.5),
		  upper_right(0.5,0.5),
		  max_radius(0.48),
		  min_radius(2e-2),
		  grid(),
		  specials_inner(),
		  specials_outer(),
		  specials_inner2(),
		  specials_inner3(),
		  specials_outer2(),
		  specials_outer3(),
		  remove_rin(),
		  remove_rout(),
	      refine_rin(),
		  refine_rout(),
		  inner_slow(),
		  outer_slow()
	  {
		  const int np=2000;
		  vector<Vector2D> MainPoints=CirclePointsRmax_a(np,min_radius,
			  max_radius,0.5*(lower_left.x+upper_right.x),0.5*(lower_left.y+
			  upper_right.y),upper_right.x,upper_right.y,lower_left.x,lower_left.y,
			  -0.25);
		  const double dr=abs(MainPoints[1]-MainPoints[0]);
		  const double dr2=abs(MainPoints[MainPoints.size()-1]-
			  MainPoints[MainPoints.size()-2]);
		  // Clear the outer points
		  vector<Vector2D> MainPoints2;
		  for(int i=0;i<(int)MainPoints.size();++i)
		  {
			  if((abs(MainPoints[i])>(min_radius+2.5*dr))&&
				  (abs(MainPoints[i])<(max_radius-2.5*dr2)))
				  MainPoints2.push_back(MainPoints[i]);
		  }
		  vector<Vector2D> inner1=Circle((int)(2*M_PI*min_radius/dr+1),min_radius+1.5*dr);
		  vector<Vector2D> inner2=Circle((int)(2*M_PI*min_radius/dr+1),min_radius+0.15*dr);
		  vector<Vector2D> inner3=Circle((int)(2*M_PI*min_radius/dr+1),min_radius);
		  vector<Vector2D> outer1=Circle((int)(2*M_PI*max_radius/dr2),max_radius-1.5*dr2);
		  vector<Vector2D> outer2=Circle((int)(2*M_PI*max_radius/dr2),max_radius-0.15*dr2);
		  vector<Vector2D> outer3=Circle((int)(2*M_PI*max_radius/dr2),max_radius);
		  specials_inner3 = (int)inner3.size();
		  specials_inner2 = (int)inner2.size();
		  specials_inner = (int)inner1.size();
		  specials_outer3 = (int)outer3.size();
		  specials_outer2 = (int)outer2.size();
		  specials_outer = (int)outer1.size();
		  grid=join(inner3,join(outer3,join(inner2,join(outer2,
			  join(inner1,join(outer1,MainPoints2))))));
		  remove_rin=min_radius+2*dr;
		  remove_rout=max_radius-2*dr2;
		  refine_rin=min_radius+4*dr;
		  refine_rout=max_radius-3.5*dr2;
		  inner_slow=min_radius+dr;
		  outer_slow=max_radius-dr2;
	  }

	  const Vector2D lower_left;
	  const Vector2D upper_right;
	  const double max_radius;
	  const double min_radius;
	  vector<Vector2D> grid;
	  int specials_inner;
	  int specials_outer;
	  int specials_inner2;
	  int specials_inner3;
	  int specials_outer2;
	  int specials_outer3;
	  double remove_rin;
	  double remove_rout;
	  double refine_rin;
	  double refine_rout;
	  double inner_slow;
	  double outer_slow;
};

namespace {
	void my_main_loop(hdsim& sim)
	{
		const double endtime=10;
		ConsecutiveSnapshots diag(endtime/1000);
		try{
			double MaxVol=1.5/sim.GetCellNo();
			GridData gd;
			const int total_special = gd.specials_inner+gd.specials_inner2+
				gd.specials_outer+gd.specials_outer2+gd.specials_inner3+
				gd.specials_outer3;
			/*			const int outer_loc = gd.specials_inner3+gd.specials_inner2+
						gd.specials_outer3+gd.specials_inner;*/
			DiskRemove disk_remove
				(gd.remove_rin,
				gd.remove_rout,
				total_special);
			DiskRefine disk_refine
				(total_special,
				1e-3,
				gd.refine_rin,
				gd.refine_rout,
				MaxVol);
			sim.SetEndTime(endtime);
			while(sim.GetTime()<endtime){
				sim.TimeAdvance2Mid();

				const vector<int> removed = sim.RemoveCells(&disk_remove);
				sim.RefineCells(&disk_refine,removed);

				diag.diagnose(sim);
			}
		}
		catch(UniversalError& eo){
			eo.AddEntry("time",sim.GetTime());
			throw;
		}
	}
}

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

	SimData(void): gd(),mm(),tess_(), eos_(mm.adiabatic_index),pm_naive_(),
		pm_inter_(pm_naive_,0.75,0.05),point_motion_(pm_inter_,gd.inner_slow,
		gd.outer_slow,gd.specials_inner+gd.specials_inner2+gd.specials_inner3+
		gd.specials_outer+gd.specials_outer2+gd.specials_outer3),eulerian_(),rs_(),
		acceptor_gravity_acceleration_(mm.acceptor_mass,0,mm.center),
		acceptor_gravity_(&acceptor_gravity_acceleration_),
		donor_gravity_acceleration_(mm.donor_mass,0,mm.donor_center),
		donor_gravity_(&donor_gravity_acceleration_),centripetal_acceleration_
		(mm.angular_velocity,mm.center_of_gravity),centripetal_force_
		(&centripetal_acceleration_),coriolis_(mm.angular_velocity),
		force_(four_vector<SourceTerm*>(&donor_gravity_,&acceptor_gravity_,
		&centripetal_force_,&coriolis_)),outer_(gd.lower_left.x,
		gd.upper_right.y,gd.upper_right.x,gd.lower_left.y),hbc_(rs_),naive_interp_
		(outer_,&hbc_,true,false),interp_(naive_interp_,eos_),sinterp_(&interp_,
		&hbc_),
		sim_(gd.grid,&tess_,&interp_,Uniform2D(mm.density),Uniform2D(mm.pressure),
		Uniform2D(0),Uniform2D(0),eos_,rs_,&point_motion_,&force_,&outer_,&hbc_,
		false,false),
		outflow_(Ratchet::in,sim_.GetAllCells()),inflow_(CalcPrimitive(
		mm.wind_density,mm.wind_pressure,mm.wind_velocity,eos_))
	  {
		  pm_inter_.SetColdFlows(&hbc_);
		  sim_.SetColdFlows(0.05,0.05);
		  sim_.setTracerInterpolation(&sinterp_);
		  sim_.SetCfl(0.8);

		  for(int i=0;i<gd.specials_inner;++i)
			  sim_.CellsEvolve[i]=&outflow_;
		  for(int i=gd.specials_inner;
			  i<gd.specials_inner+gd.specials_outer;++i){
				  const Vector2D pos = sim_.GetMeshPoint(i);
				  if(abs(pos.y)<mm.orifice&&pos.x>0)
					  sim_.CellsEvolve[i] = &inflow_;
				  else
					  sim_.CellsEvolve[i] = &outflow_;
		  }
	  }

	  hdsim& getSim(void)
	  {
		  return sim_;
	  }

	  GridData gd;
	  MagicNumbers mm;

private:
	VoronoiMesh tess_;
	const IdealGas eos_;
	Lagrangian pm_naive_;
	RoundCells pm_inter_;
	EdgeRepellant point_motion_;
	Eulerian eulerian_;
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
	RigidWallHydro hbc_;
	LinearGauss naive_interp_;
	EOSConsistent interp_;
	LinearGaussScalar sinterp_;
	hdsim sim_;
	Ratchet outflow_;
	ConstantPrimitive inflow_;
};

class SimDataReset
{
public:
	SimDataReset(ResetDump const& dump): gd(),mm(),tess_(), eos_(mm.adiabatic_index),
		pm_naive_(),pm_inter_(pm_naive_,0.75,0.05),
		point_motion_(pm_inter_,gd.inner_slow,
		gd.outer_slow,gd.specials_inner+gd.specials_inner2+gd.specials_inner3+
		gd.specials_outer+gd.specials_outer2+gd.specials_outer3),eulerian_(),rs_(),
		acceptor_gravity_acceleration_(mm.acceptor_mass,0,mm.center),
		acceptor_gravity_(&acceptor_gravity_acceleration_),
		donor_gravity_acceleration_(mm.donor_mass,0,mm.donor_center),
		donor_gravity_(&donor_gravity_acceleration_),centripetal_acceleration_
		(mm.angular_velocity,mm.center_of_gravity),centripetal_force_
		(&centripetal_acceleration_),coriolis_(mm.angular_velocity),
		force_(four_vector<SourceTerm*>(&donor_gravity_,&acceptor_gravity_,
		&centripetal_force_,&coriolis_)),outer_(gd.lower_left.x,
		gd.upper_right.y,gd.upper_right.x,gd.lower_left.y),hbc_(rs_),naive_interp_
		(outer_,&hbc_,true,false),interp_(naive_interp_,eos_),sinterp_(&interp_,
		&hbc_),
		sim_(dump.snapshot.mesh_points,&tess_,&interp_,dump.snapshot.cells,
		eos_,rs_,&point_motion_,&force_,&outer_,&hbc_,dump.tracers,dump.time,
		dump.cfl,dump.cycle,dump.coldflows,dump.a,dump.b,dump.densityfloor,
		dump.densitymin,dump.pressuremin),
		outflow_(Ratchet::in,sim_.GetAllCells()),inflow_(CalcPrimitive(
		mm.wind_density,mm.wind_pressure,mm.wind_velocity,eos_))
	  {
		  pm_inter_.SetColdFlows(&hbc_);
		  sim_.setTracerInterpolation(&sinterp_);

		  for(int i=0;i<gd.specials_inner3;++i)
			  sim_.CellsEvolve[i]=&outflow_;
		  for(int i=gd.specials_inner3;
			  i<gd.specials_inner3+gd.specials_outer3;++i){
				  const Vector2D pos = sim_.GetMeshPoint(i);
				  if(abs(pos.y)<mm.orifice&&pos.x>0)
					  sim_.CellsEvolve[i] = &inflow_;
				  else
					  sim_.CellsEvolve[i] = &outflow_;
		  }
	  }

	  hdsim& getSim(void)
	  {
		  return sim_;
	  }

	  GridData gd;
	  MagicNumbers mm;

private:
	VoronoiMesh tess_;
	const IdealGas eos_;
	Lagrangian pm_naive_;
	RoundCells pm_inter_;
	EdgeRepellant point_motion_;
	Eulerian eulerian_;
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
	RigidWallHydro hbc_;
	LinearGauss naive_interp_;
	EOSConsistent interp_;
	LinearGaussScalar sinterp_;
	hdsim sim_;
	Ratchet outflow_;
	ConstantPrimitive inflow_;
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
	/*const IdealGas eos(1.2);
	ResetDump dump;
	read_hdf5_snapshot(dump,"snapshot_33.h5",&eos);
	ResetRead("c:\\reset.bin",dump,&eos);
	SimDataReset sim_data(dump);*/
  
	SimData sim_data;
	hdsim& sim = sim_data.getSim();
	sim.SetCfl(0.8);

	write_snapshot_to_hdf5(sim,"initial.h5"); 
	//ResetOutput("c:\\reset.bin",sim);
	try{
		my_main_loop(sim);
	}
	catch(UniversalError const& eo){
		report_error(eo,sim);
		throw;
	}

	write_snapshot_to_hdf5(sim,"final.h5");
}
