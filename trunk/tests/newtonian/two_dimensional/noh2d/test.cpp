#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/CustomOuter.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/misc/int2str.hpp"
#include "source/newtonian/test_2d/noh2d/noh_hbc.hpp"
#include "source/newtonian/test_2d/noh2d/noh_amr.hpp"

class xVel :public SpatialDistribution
{
private:
	double v_;
public:
	xVel(double v):v_(v){};
	~xVel(){};
	double EvalAt(Vector2D const& point) const
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
	double EvalAt(Vector2D const& point) const
	{
		double r=abs(point);
		return -v_*point.y/r;
	}
};


int main(void)
{
	// Set up the initial grid points
	int np = read_int("resolution.txt");
	double width=2;
	vector<Vector2D> InitPoints=SquareMesh(np,np,width,width);

	// Set up the boundary type for the points
	SquareBox outer(-width*0.5,width*0.5,width*0.5,-width*0.5);

	// Set up the tessellation
	VoronoiMesh tess;

	// Set up the Riemann solver
	Hllc rs;

	// Set the hydro boundary conditions
	Vector2D center(0,0);
	double rho=1;
	double v_in=1;
	double p=1e-6;
	NohHBC hbc(center,rho,v_in,p,rs);

	// Set up the equation of state
	IdealGas eos(read_number("adiabatic_index.txt"));

	// Set up the point motion scheme
	Lagrangian l_motion;
	RoundCells pointmotion(l_motion,hbc);
	pointmotion.SetColdFlows();

	// Set up the interpolation
	LinearGaussConsistent interpolation(eos,outer,&hbc);

	// Set up the initial Hydro
	Uniform2D density(rho);
	Uniform2D pressure(p);
	xVel xvelocity(v_in);
	yVel yvelocity(v_in);

	// Set up the external source term
	ZeroForce force;

	// Set up the simulation
	hdsim sim(InitPoints,&tess,&interpolation,density,pressure,xvelocity,
		yvelocity,eos,rs,&pointmotion,&force,&outer,&hbc);

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
	double cfl=0.3;
	sim.SetCfl(cfl);

	// How long shall we run the simulation?
	double tend=2;
	sim.SetEndTime(tend);

	// Run main loop of the sim
	while(sim.GetTime()<tend)
	{
		try
		{
			// Advance one time step
			sim.TimeAdvance2Mid();
			// Remove small cells
			vector<int> removed_cells=sim.RemoveCells(&remove);
			// Refine big cells
			sim.RefineCells(&refine,removed_cells);
		}
		catch(UniversalError const& eo)
		{
			DisplayError(eo,sim.GetCycle());
		}
	}

	// Done running the simulation, output the data
	string filename="final.h5";
	write_snapshot_to_hdf5(sim,filename);

	return 0;
}
