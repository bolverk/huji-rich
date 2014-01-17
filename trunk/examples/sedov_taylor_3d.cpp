#include "../source/tessellation/VoronoiMesh.hpp"
#include "../source/newtonian/two_dimensional/hdsim2d.hpp"
#include "../source/newtonian/common/hllc.hpp"
#include "../source/newtonian/common/ideal_gas.hpp"
#include "../source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "../source/newtonian/two_dimensional/spatial_distributions/Circle2D.hpp"
#include "../source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "../source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "../source/newtonian/two_dimensional/source_terms/cylinderical_geometry.hpp"
#include "../source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "../source/newtonian/two_dimensional/linear_gauss_scalar.hpp"
#include "../source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "../source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "../source/misc/mesh_generator.hpp"
#include "../source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "../source/misc/int2str.hpp"

int main(void)
{
	// Set up the initial grid points
	int npointsx=50;
	int npointsy=50;
	double widthx=1;
	double widthy=1;
	vector<Vector2D> InitPoints=SquareMesh(npointsx,npointsy,widthx,widthy,false);

	// Set up the boundary type for the points
	SquareBox outer(0,widthx,widthy,0);

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
	double rho=1;
	double low_pressure=1;
	double high_pressure=400;
	double high_radius=0.15;
	double x_velocity=0;
	double y_velocity=0;
	Uniform2D density(rho);
	Circle2D pressure(0,0,high_radius,high_pressure,low_pressure);
	Uniform2D xvelocity(x_velocity);
	Uniform2D yvelocity(y_velocity);

	// Set up the external source term
	CylindericalGeometry force(Vector2D(0,0),Vector2D(0,1));

	// Set up the simulation
	hdsim sim(InitPoints,&tess,&interpolation,density,pressure,xvelocity,
		yvelocity,eos,rs,&pointmotion,&force,&outer,&hbc);

	// Add the tracer for following the mass
	Circle2D tracer(0,0,high_radius,1,0);
	sim.addTracer(tracer);

	// Add tracer interpolation
	LinearGaussScalar sinterp(&interpolation,&hbc);
	sim.setTracerInterpolation(&sinterp);

	// Choose the Courant number
	double cfl=0.7;
	sim.SetCfl(cfl);

	// How long shall we run the simulation?
	double tend=0.1;
	sim.SetEndTime(tend);

	// Custom output criteria
	double output_dt=0.01;
	double last_dump_time=0;
	int dump_number=0;
	write_snapshot_to_hdf5(sim,"c:\\sim_data\\output"+
		int2str(dump_number)+".bin");
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
	write_snapshot_to_hdf5(sim,filename);
	
	// We are done!!
	cout<<"Finished running the simulation"<<endl;

	return 0;
}