#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/PeriodicBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/PeriodicHydro.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
//#include "source/newtonian/two_dimensional/interpolations/LinearGaussArepo.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/misc/int2str.hpp"
#include "source/newtonian/two_dimensional/RefineStrategy.hpp"
#include "source/newtonian/two_dimensional/RemovalStrategy.hpp"
#include "source/tessellation/voronoi_logger.hpp"

class LeftRemove: public RemovalStrategy
{
public:
	LeftRemove(double xl,double minA):xl_(xl),minA_(minA){}
	vector<int> CellsToRemove(Tessellation const& tess,
		vector<Primitive> const& /*cells*/,vector<vector<double> > const& /*tracers*/,
		double /*time*/)const
	{
		int n=tess.GetPointNo();
		vector<int> res;
		vector<double> merits;
		for(int i=0;i<n;++i)
		{
			const double v=tess.GetVolume(i);
			const Vector2D point=tess.GetMeshPoint(i);
			if(point.x<xl_&&v<minA_)
			{
				res.push_back(i);
				merits.push_back(-point.x);
			}
		}
		res=RemoveNeighbors(merits,res,tess);
		CheckOutput(tess,res);
		return res;
	}
private:
	const double xl_,minA_;
};

class RightAdd: public RefineStrategy
{
public:
	RightAdd(double xr,double maxA):xr_(xr),maxA_(maxA){}
	vector<int> CellsToRefine(Tessellation const& tess,
		vector<Primitive> const& /*cells*/,vector<vector<double> > const& /*tracers*/,
		double /*time*/,vector<Vector2D> &directions,vector<int> const& Removed)
	{
		vector<int> res;
		int n=tess.GetPointNo();
		directions.clear();
		Vector2D dir(1,0);
		for(int i=0;i<n;++i)
		{
			const double v=tess.GetVolume(i);
			const Vector2D point=tess.GetMeshPoint(i);
			if(point.x>xr_&&v>maxA_)
			{
				res.push_back(i);
				//directions.push_back(dir);
			}
		}
		return RemoveDuplicatedLately(res,n,directions,Removed,tess);
	}
private:
	double xr_,maxA_;

};

namespace 
{

	double GetMass(hdsim const& sim)
	{
		int n=sim.GetCellNo();
		double res=0;
		for(int i=0;i<n;++i)
			res+=sim.GetCellVolume(i)*sim.GetCell(i).Density;
		return res;
	}

	Vector2D GetMomentum(hdsim const& sim)
	{
		int n=sim.GetCellNo();
		Vector2D res;
		for(int i=0;i<n;++i)
			res+=sim.GetCellVolume(i)*sim.GetCell(i).Density*sim.GetCell(i).Velocity;
		return res;
	}

	double GetMaxDensity(hdsim const& sim)
	{
		int n=sim.GetCellNo();
		double res=0;
		for(int i=0;i<n;++i)
			res=max(res,sim.GetCell(i).Density);
		return res;
	}

	double GetMinDensity(hdsim const& sim)
	{
		int n=sim.GetCellNo();
		double res=sim.GetCell(0).Density;
		for(int i=1;i<n;++i)
			res=min(res,sim.GetCell(i).Density);
		return res;
	}

}

int main(void)
{
	// Set up the initial grid points
	int np = 30;
	double width=1;
	vector<Vector2D> InitPoints=SquareMesh(np,np,width,width);

	// Set up the boundary type for the points
	PeriodicBox outer(-width*0.5,width*0.5,width*0.5,-width*0.5);

	// Set up the tessellation
	VoronoiMesh tess;

	// Set up the Riemann solver
	Hllc rs;

	// Set the hydro boundary conditions
	PeriodicHydro hbc(rs);

	// Set up the equation of state
	IdealGas eos(5./3.);

	// Set up the point motion scheme
	Lagrangian l_motion;
	RoundCells pointmotion(l_motion,hbc,0.6);

	// Set up the interpolation
	LinearGaussConsistent interpolation(eos,outer,hbc);

	// Set up the initial Hydro
	Uniform2D density(1);
	Uniform2D pressure(1);
	Uniform2D xvelocity(1);
	Uniform2D yvelocity(0);

	// Set up the external source term
	ZeroForce force;

	// Set up the simulation
	hdsim sim(InitPoints,tess,interpolation,density,pressure,xvelocity,
		yvelocity,eos,rs,pointmotion,force,outer,hbc);

	// Define the AMR classes
	double Vmax=0.8*width*width/(np*np);
	double Vmin=1.2*width*width/(np*np);
	LeftRemove remove(-width*0.4,Vmin);
	RightAdd refine(width*0.4,Vmax);

	// How long shall we run the simulation?
	double tend=1.5;
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
			DisplayError(eo);
		}
	}

	// Run diagnostics
	vector<double> result;
	result.push_back(GetMass(sim));
	result.push_back(GetMomentum(sim).x);
	result.push_back(GetMomentum(sim).y);
	result.push_back(GetMinDensity(sim));
	result.push_back(GetMaxDensity(sim));

	write_vector(result,"result.txt");
	return 0;
}
