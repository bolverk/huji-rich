#include "../source/tessellation/VoronoiMesh.hpp"
#include "../source/newtonian/two_dimensional/hdsim2d.hpp"
#include "../source/newtonian/common/hllc.hpp"
#include "../source/newtonian/common/ideal_gas.hpp"
#include "../source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "../source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "../source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "../source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "../source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "../source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "../source/misc/mesh_generator.hpp"
#include "../source/newtonian/two_dimensional/hdf5_diagnostics.hpp"

int main(void)
{
	// Set up the initial grid size
	double widthx=2;
	double widthy=2;

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

	// Set up the external source term
	ZeroForce force;

	// Read the snapshot file
	ResetDump dump;
	read_hdf5_snapshot(dump,"c:\\sim_data\\output.bin",&eos);

	// Set up the simulation
	hdsim sim(dump.snapshot.mesh_points,&tess,&interpolation,dump.snapshot.cells,
		eos,rs,&pointmotion,&force,&outer,&hbc,dump.tracers,dump.time,dump.cfl,
		dump.cycle,dump.coldflows,dump.a,dump.b,dump.densityfloor,dump.densitymin,
		dump.pressuremin);

	// How long shall we run the simulation?
	double tend=0.1;
	sim.SetEndTime(tend);

	// Run main loop of the sim
	while(sim.GetTime()<tend)
	{
		// This purely for user feedback
		if(sim.GetCycle()%25==0)
			cout<<"Sim time is "<<sim.GetTime()<<" Step number "<<sim.GetCycle()<<endl;
		// Advance one time step
		sim.TimeAdvance2Mid();
	}

	// Done running the simulation, output the data
	string filename="c:\\sim_data\\output.bin";
	write_snapshot_to_hdf5(sim,filename);
	
	// We are done!!
	cout<<"Finished running the simulation"<<endl;

	return 0;
}