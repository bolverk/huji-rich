#include "../source/tessellation/VoronoiMesh.hpp"
#include "../source/newtonian/two_dimensional/hdsim2d.hpp"
#include "../source/newtonian/common/hllc.hpp"
#include "../source/newtonian/common/ideal_gas.hpp"
#include "../source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "../source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "../source/newtonian/two_dimensional/hydro_boundary_conditions/InFlow.hpp"
#include "../source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "../source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "../source/misc/mesh_generator.hpp"
#include "../source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "../source/misc/int2str.hpp"
#include "../source/newtonian/two_dimensional/source_terms/CenterGravity.hpp"
#include "../source/newtonian/two_dimensional/source_terms/SeveralSources.hpp"
#include "../source/newtonian/two_dimensional/custom_evolutions/ConstantPrimitiveEvolution.hpp"
#include "../source/newtonian/two_dimensional/linear_gauss_scalar.hpp"
#include "binary_outflow.hpp"

vector<int> FindCells(Vector2D const& point,Tessellation const& tess,int nx)
{
	int N=tess.GetPointNo();
	int cell_index=0;
	double distance=tess.GetMeshPoint(0).distance(point);
	for(int i=1;i<N;++i)
	{
		if(point.distance(tess.GetMeshPoint(i))<distance)
		{
			distance=point.distance(tess.GetMeshPoint(i));
			cell_index=i;
		}
	}
	vector<int> result;
	result.push_back(cell_index);
	result.push_back(cell_index+1);
	result.push_back(cell_index-1);
	result.push_back(cell_index+nx);
	result.push_back(cell_index-nx);
	result.push_back(cell_index+nx-1);
	result.push_back(cell_index-nx-1);
	result.push_back(cell_index+nx+1);
	result.push_back(cell_index-nx+1);
	return result;
}

int main(void)
{
	// Set up the initial grid points
	int npointsx=100;
	int npointsy=100;
	double widthx=2;
	double widthy=2;
	vector<Vector2D> InitPoints=SquareMesh(npointsx,npointsy,widthx,widthy);

	// Set up the boundary type for the points
	SquareBox outer(-widthx/2,widthx/2,widthy/2,-widthy/2);

	// Set up the tessellation
	VoronoiMesh tess;

	// Set up the Riemann solver
	Hllc rs;

	// Set up the equation of state
	double gamma=5./3.;
	IdealGas eos(gamma);

	// Set up the point motion scheme
	Eulerian pointmotion;

	// Set up the initial Hydro
	double rho=0.001;
	double P=0.001;
	Uniform2D density(rho);
	Uniform2D pressure(P);
	Uniform2D xvelocity(0);
	Uniform2D yvelocity(0);

	
	// Set the hydro boundary conditions
	Primitive inflow_primitive;
	inflow_primitive.Density=rho;
	inflow_primitive.Pressure=P;
	inflow_primitive.Velocity=Vector2D(0,0);
	inflow_primitive.Energy=eos.dp2e(rho,P);
	inflow_primitive.SoundSpeed=eos.dp2c(rho,P);
	vector<double> outer_entropy;
	outer_entropy.push_back(eos.dp2s(rho,P));
	InFlow hbc(inflow_primitive,rs,outer_entropy);

	// Set up the interpolation
	LinearGaussConsistent interpolation(eos,outer,&hbc);

	// Set up the external source term
	double star_mass=1;
	double softening_scale=0.01;
	double x1=-0.4;
	double x2=0.3;
	double y1=0;
	double y2=0;
	double mdot=0.01;
	double cs=9.5;
	CenterGravity gravity1(star_mass,softening_scale,Vector2D(x1,y1));
	CenterGravity gravity2(4*star_mass,softening_scale,Vector2D(x2,y2));
	BinaryOutflow star_output(mdot,0.01*gamma*cs*cs/(gamma-1),cs,npointsx,npointsy);
	ConservativeForce gforce1(&gravity1);
	ConservativeForce gforce2(&gravity2);
	vector<SourceTerm*> forces;
	forces.push_back(&gforce1);
	forces.push_back(&gforce2);
	forces.push_back(&star_output);

	SeveralSources force(forces);

	// Set up the simulation
	hdsim sim(InitPoints,&tess,&interpolation,density,pressure,xvelocity,
		yvelocity,eos,rs,&pointmotion,&force,&outer,&hbc);

	// Set cold flows on
	double kineticfraction=0.01;
	double gravityfraction=0.01;
	sim.SetColdFlows(kineticfraction,gravityfraction);

	// Define the interpolation for the cold flows
	LinearGaussScalar sinterp(&interpolation,&hbc);
	sim.setTracerInterpolation(&sinterp);

	// Initialize the outflow
	star_output.Init(tess,x1,y1);

	// Set the outflow cell to be constant
	ConstantPrimitiveEvolution ouflow_evolve;
	sim.CellsEvolve[star_output.GetCellIndex()]=&ouflow_evolve;

	// Set custom evolution for 2nd star in order to remove mass
	vector<int> sink_cells=FindCells(Vector2D(x2,y2),tess,npointsx);
	for(size_t i=0;i<sink_cells.size();++i)
		sim.CellsEvolve[sink_cells[i]]=&ouflow_evolve;

	// Set the first time step
	sim.SetTimeStepExternal(0.001);

	// Choose the Courant number
	double cfl=0.7;
	sim.SetCfl(cfl);

	// How long shall we run the simulation?
	double tend=1;
	sim.SetEndTime(tend);

	// Custom output criteria
	double output_dt=0.01;
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
				BinOutput("c:\\sim_data\\output"+int2str(dump_number)+".bin",sim,tess);
			}

			// Advance one time step
			double external_dt=min(gforce1.GetTimeStep(),gforce2.GetTimeStep());
			if(sim.GetCycle()>0)
				sim.SetTimeStepExternal(external_dt);
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