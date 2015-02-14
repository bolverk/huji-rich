#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/two_dimensional/Reset.hpp"
#include "source/misc/int2str.hpp"
#include "source/tessellation/RoundGrid.hpp"
#include "source/newtonian/test_2d/noh2d/noh_amr.hpp"
#include "source/newtonian/test_2d/noh2d/noh_hbc.hpp"

class xVel :public SpatialDistribution
{
private:
	double v_;
public:
	xVel(double v):v_(v){};
	~xVel(){};
	double operator()(Vector2D const& point) const
	{
		double r=abs(point);
		return -v_*point.x/r;;
	}
};

class yVel :public SpatialDistribution
{
private:
	double v_;
public:
	yVel(double v):v_(v){};
	~yVel(){};
	double operator()(Vector2D const& point) const
	{
		double r=abs(point);
		return -v_*point.y/r;
	}
};


int main(void)
{
	// Set up the initial grid points
	int np = 50;
	double width=1;
	vector<Vector2D> InitPoints=RandSquare(np*np,-width/2,width/2,-width/2,width/2);

	// Set up the boundary type for the points
	SquareBox outer(-width/2,width/2,width/2,-width/2);

	// Set up the tessellation
	VoronoiMesh tess;

	// Set up the Riemann solver
	Hllc rs;

	// Set the hydro boundary conditions
	Vector2D center(0,0);
	double rho=1;
	double v_in=1;
	double p=1e-6;
	NohHBC hbc(center,rho,v_in,p);

	// Set up the equation of state
	double gamma=5./3.;
	IdealGas eos(gamma);

	// Set up the point motion scheme
	Lagrangian l_motion;
	RoundCells pointmotion(l_motion,hbc,0.7,0.02,true);

	// Set up the interpolation
	LinearGaussConsistent interpolation(eos,outer,hbc);

	// Set up the initial Hydro
	Uniform2D density(rho);
	Uniform2D pressure(p);
	xVel xvelocity(v_in);
	yVel yvelocity(v_in);

	// Set up the external source term
	ZeroForce force;

	// Make the grid nice and round for first time step
	InitPoints=RoundGrid(InitPoints,&outer,5);

	// Set up the simulation
	hdsim sim(InitPoints,tess,interpolation,density,pressure,xvelocity,
		yvelocity,eos,rs,pointmotion,force,outer,hbc);

	// Set cold flows on
	double kineticfraction=0.01;
	double gravityfraction=0.01;
	sim.SetColdFlows(kineticfraction,gravityfraction);

	// Define the AMR classes
	double Vmax=2*width*width/(np*np);
	double Vmin=0.25*width*width/(np*np);
	NohRefine refine(Vmax);
	NohRemove remove(Vmin);
	
	// Choose the Courant number
	double cfl=0.7;
	sim.SetCfl(cfl);

	// How long shall we run the simulation?
	double tend=0.5;
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
			// Remove small cells
			vector<int> removed_cells=sim.RemoveCells(&remove);
			// Refine big cells
			sim.RefineCells(&refine,removed_cells);
		}
		catch(UniversalError const& eo)
		{
			DisplayError(eo);
		}
	}

	// Done running the simulation, output the data
	string filename="c:\\sim_data\\output.bin";
	write_snapshot_to_hdf5(sim,filename);

	// We are done!!
	cout<<"Finished running the simulation"<<endl;


	return 0;
}
