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

using namespace std;
using namespace simulation2d;

class MagicNumbers
{
public:

	MagicNumbers(void):
	  center(0,0),
		  adiabatic_index(5./3.),
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
		  orifice(0.01),
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
		  min_radius(1e-2),
		  gsq(1.005),
		  grid(),
		  specials_inner(),
		  specials_outer(),
		  specials_inner2(),
		specials_inner3(),
		  specials_outer2(),
		  ring_radii()
	  {
		  const double rr = max_radius/min_radius;
		  vector<double> radii_list(int(log(1+(gsq-1)*rr)/log(gsq)));
		  for(int i=0;i<(int)radii_list.size();++i)
			  radii_list[i] = min_radius*(pow(gsq,i+1)-1)/(gsq-1);
		  vector<vector<Vector2D> > rings(radii_list.size(),vector<Vector2D>());
		  {
			  const int nq = int(2*M_PI);
			  for(int j=0;j<nq;++j)
				  rings[0].push_back(polar2cart(radii_list[0],
				  2*M_PI*double(j)/double(nq)));
		  }
		  for(int i=1;i<(int) radii_list.size();++i){
			  const int nq = int(2*M_PI*radii_list[i]/
				  (radii_list[i]-radii_list[i-1]));
			  for(int j=0;j<nq;++j)
				  rings[i].push_back(polar2cart(radii_list[i],
				  2*M_PI*double(j)/double(nq)));
		  }

		  specials_inner = (int)rings[0].size();
		  specials_outer = (int)rings[rings.size()-1].size();
		  specials_inner2 = (int)rings[1].size();
		  specials_inner3 = (int)rings[2].size();
		  specials_outer2 = (int)rings[rings.size()-2].size();
		  for(int j=0;j<specials_inner;++j)
			  grid.push_back(rings[0][j]);
		  for(int j=0;j<specials_outer;++j)
			  grid.push_back(rings[rings.size()-1][j]);
		  for(int j=0;j<specials_inner2;++j)
			  grid.push_back(rings[1][j]);
		  for(int j=0;j<specials_inner3;++j)
			  grid.push_back(rings[2][j]);
		  for(int j=0;j<specials_outer2;++j)
			  grid.push_back(rings[rings.size()-2][j]);
		  for(int i=3;i<(int)rings.size()-2;++i){
			  for(int j=0;j<(int)rings[i].size();++j)
				  grid.push_back(rings[i][j]);
		  } 
		  ring_radii = radii_list;
	  }

	  const Vector2D lower_left;
	  const Vector2D upper_right;
	  const double max_radius;
	  const double min_radius;
	  const double gsq;
	  vector<Vector2D> grid;
	  int specials_inner;
	  int specials_outer;
	  int specials_inner2;
	  int specials_inner3;
	  int specials_outer2;
	  vector<double> ring_radii;
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
				gd.specials_outer+gd.specials_outer2+gd.specials_inner3;
			DiskRemove disk_remove
				(gd.ring_radii[3],
				gd.ring_radii[gd.ring_radii.size()-3],
				total_special);
			DiskRefine disk_refine
				(total_special,
				1e-3,
				gd.ring_radii[6],
				MaxVol);
			cout<<"Refine ring "<<gd.ring_radii[6]<<endl;
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

	SimData(void): gd(),mm(),tess_(), eos_(5./3.),pm_naive_(),
		pm_inter_(pm_naive_,0.75,0.05),point_motion_(pm_inter_,gd.ring_radii[3],
		gd.ring_radii[gd.ring_radii.size()-3],gd.specials_inner,gd.specials_outer,
		gd.specials_inner2,gd.specials_inner3,gd.specials_outer2,
		sqrt(mm.acceptor_mass/pow(gd.ring_radii[0],3))),eulerian_(),rs_(),
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
	SimDataReset(ResetDump const& dump): gd(),mm(),tess_(), eos_(5./3.),pm_naive_(),
		pm_inter_(pm_naive_,0.75,0.05),point_motion_(pm_inter_,gd.ring_radii[3],
		gd.ring_radii[gd.ring_radii.size()-3],gd.specials_inner,gd.specials_outer,
		gd.specials_inner2,gd.specials_inner3,gd.specials_outer2,
		sqrt(mm.acceptor_mass/pow(gd.ring_radii[0],3))),eulerian_(),rs_(),
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
	/*const IdealGas eos(5./3.);
	ResetDump dump;
	read_hdf5_snapshot(dump,"snapshot_103.h5",&eos);
	SimDataReset sim_data(dump);
	*/
	SimData sim_data;
	
	hdsim& sim = sim_data.getSim();

	write_snapshot_to_hdf5(sim,"initial.h5"); 

	try{
		my_main_loop(sim);
	}
	catch(UniversalError const& eo){
		report_error(eo,sim);
		throw;
	}

	write_snapshot_to_hdf5(sim,"final.h5");
}
