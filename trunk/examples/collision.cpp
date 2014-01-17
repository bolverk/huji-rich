#include "../source/tessellation/VoronoiMesh.hpp"
#include "../source/newtonian/two_dimensional/hdsim2d.hpp"
#include "../source/newtonian/common/hllc.hpp"
#include "../source/newtonian/common/ideal_gas.hpp"
#include "../source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "../source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "../source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "../source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "../source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "../source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "../source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "../source/misc/mesh_generator.hpp"
#include "../source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "../source/misc/int2str.hpp"
#include "collide_density.hpp"

int main(void)
{
	// Set up the initial grid points
	int npointsx=50;
	int npointsy=50;
	double widthx=2;
	double widthy=2;
	double blob_radius=0.2;
	int npoints_blob=1000;
	double blob1_xc=-0.5;
	double blob2_xc=0.5;
	double blob_yc=0;

	vector<Vector2D> background=SquareMesh(npointsx,npointsy,widthx,widthy);
	vector<Vector2D> blob1=CirclePointsRmax(npoints_blob,blob_radius*0.1,
		blob_radius,blob1_xc,blob_yc);
	vector<Vector2D> blob2=CirclePointsRmax(npoints_blob,blob_radius*0.1,
		blob_radius,blob2_xc,blob_yc);
	vector<Vector2D> InitPoints=background;
	InitPoints.insert(InitPoints.end(),blob1.begin(),blob1.end());
	InitPoints.insert(InitPoints.end(),blob2.begin(),blob2.end());

	// Set up the boundary type for the points
	SquareBox outer(-widthx/2,widthx/2,widthy/2,-widthy/2);

	// Set up the tessellation
	VoronoiMesh tess;

	// Set up the Riemann solver
	Hllc rs;

	// Set the hydro boundary conditions
	RigidWallHydro hbc(rs);

	// Set up the equation of state
	double gamma=5./3.;
	IdealGas eos(gamma);

	// Set up the point motion scheme
	Lagrangian l_motion;
	RoundCells pointmotion(l_motion);

	// Set up the interpolation
	LinearGaussConsistent interpolation(eos,outer,&hbc);

	// Set up the initial Hydro
	double density1=10;
	double density2=10;
	double background_density=0.1;
	double fluid_pressure=1;
	double xvelocity1=1;
	double xvelocity2=-1;
	double background_velocity=0;

	CollideDensity density(density1,density2,background_density,blob1_xc,blob_yc,
		blob2_xc,blob_yc,blob_radius,blob_radius);
	Uniform2D pressure(fluid_pressure);
	CollideDensity xvelocity(xvelocity1,xvelocity2,background_velocity,blob1_xc,
		blob_yc,blob2_xc,blob_yc,blob_radius,blob_radius);
	Uniform2D yvelocity(background_velocity);

	// Set up the external source term
	ZeroForce force;

	// Set up the simulation
	hdsim sim(InitPoints,&tess,&interpolation,density,pressure,xvelocity,
		yvelocity,eos,rs,&pointmotion,&force,&outer,&hbc);

	// Choose the Courant number
	double cfl=0.7;
	sim.SetCfl(cfl);

	// How long shall we run the simulation?
	double tend=1;
	sim.SetEndTime(tend);

	// Custom output criteria
	double output_dt=0.1;
	double last_dump_time=0;
	int dump_number=0;

	// Run main loop of the sim
	while(sim.GetTime()<tend)
	{
		try
		{
			// This purely for user feedback
			if(sim.GetCycle()%25==0)
				cout<<"Sim time is "<<sim.GetTime()<<" Step number "<<sim.GetCycle()<<endl;

			// Custom output criteria
			if((sim.GetTime()-last_dump_time)>output_dt)
			{
				last_dump_time=sim.GetTime();
				++dump_number;
				write_snapshot_to_hdf5(sim,"c:\\sim_data\\output"+
					int2str(dump_number)+".bin");
			}

			// Advance one time step
			sim.TimeAdvance2Mid();
		}
		catch(UniversalError const& eo)
		{
			DisplayError(eo,sim.GetCycle());
		}
	}

	// Done running the simulation, output the data
	string filename="c:\\sim_data\\output.bin";
	BinOutput(filename,sim,tess);

	// We are done!!
	cout<<"Finished running the simulation"<<endl;

	return 0;
}