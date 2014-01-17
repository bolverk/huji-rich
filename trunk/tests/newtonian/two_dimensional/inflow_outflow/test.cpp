#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include "source/misc/int2str.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/tessellation/geometry.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/tessellation/RoundGrid.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/Line2D.hpp"
#include "source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/linear_gauss_scalar.hpp"
#include "source/newtonian/common/ersig.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/CustomOuter.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/FreeFlow.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/InFlow.hpp"

#include "source/newtonian/two_dimensional/diagnostics.hpp"

using namespace std;

namespace {

void CheckSim(hdsim const& sim,Tessellation const& tess)
{
	double GoodMass=0.06;
	int n=tess.GetPointNo();
	double SimMass=0;
	for(int i=0;i<n;++i)
	{
		SimMass+=tess.GetVolume(i)*sim.GetCell(i).Density;
	}
	ofstream myfile("result.txt");
	myfile<<(GoodMass-SimMass)/GoodMass;
	myfile.close();
}

class RefineBow : public RefineStrategy
{
private:
	int _N;
public:
  RefineBow(int npoints):
    _N(npoints) {}
	~RefineBow(){};
	vector<int> CellsToRefine(Tessellation const* tess,
		vector<Primitive> const& /*cells*/,vector<vector<double> > const& /*tracers*/,
		double /*time*/,vector<Vector2D> &directions,const vector<int> &Removed)
	{
		double dv=1.2*5.0/(_N*_N);
		vector<int> res;
		directions.clear();
		Vector2D dir(1,0);
		for(int i=0;i<_N;++i)
		{
			if(tess->GetVolume(i)>dv)
			{
				res.push_back(i);
				directions.push_back(dir);
			}
		}
		return RemoveDuplicatedLately(res,tess->GetPointNo(),directions,Removed);
	}
};

class RemoveBow : public RemovalStrategy
{
public:
	vector<int> CellsToRemove(Tessellation const* tess,
		vector<Primitive> const& cells,vector<vector<double> > const& /*tracers*/,
		double /*time*/)const
	{
		vector<int> res;
		vector<double> merit;
		int n=(int)cells.size();
		double r=2/sqrt(1.0*n);
		for(int i=0;i<n;++i)
		{
			double dx=1-tess->GetMeshPoint(i).x;
			if(dx<2*r)
			{
				res.push_back(i);
				merit.push_back(1.0/dx);
			}
		}
		// Make sure there are no neighbors
		vector<int> result=RemoveNeighbors(merit,res,tess);
		CheckOutput(tess,result);
		return result;
	}
};
}

int main(void)
{
  // Initialization
  int np = 40;
  double width=1;

  SquareBox outer(-width,width,width*1.5,-width*1.5);

  vector<Vector2D> InitPoints=SquareMesh(np,np,width*2,width*3);
	
  Hllc rs;
  IdealGas eos(5./3.);
  RefineBow refine(np);
  RemoveBow remove;

  FreeFlow hbc1(rs);
  RigidWallHydro hbc3(rs);
  Primitive in;
  // Set the inflow
  in.Density=0.01;
  in.Velocity.Set(1,0);
  in.Pressure=0.0015;
  in.Energy=eos.dp2e(in.Density,in.Pressure);
  in.SoundSpeed=eos.dp2c(in.Density,in.Pressure);
  InFlow hbc2(in,rs);
  CustomOuter hbc(&hbc2,&hbc1,&hbc3,&hbc3);
  
  VoronoiMesh tess;
  LinearGaussConsistent interp(eos,outer,&hbc,true,false,0.1,0.4);

  Uniform2D pressure(0.0015);
  Uniform2D density(0.01);
  Uniform2D xvelocity(1);
  Uniform2D yvelocity(0);

  Lagrangian bpm;
  RoundCells pointmotion(bpm,hbc,0.65,0.05,np,&outer);
  pointmotion.SetColdFlows();

  ZeroForce force;

  //ERSIG rs(5./3.,"zero flux");
  hdsim sim(InitPoints, //Initial points
	    &tess, // Tessellation
	    &interp, // Spatial reconstruction
	    density, pressure, xvelocity, yvelocity,
	    eos,rs,&pointmotion,&force,&outer,&hbc);

  // Main process
  sim.SetCfl(0.8);
  int iter = 0;
  const int max_iter = 5e6;
  double tf =2;
  sim.SetEndTime(tf);
  int cntr = 0;
  string loc="c:\\sim_data\\";
  int dumpnum=1;
   while(tf>sim.GetTime())
    {
      try
	{
		if(cntr%20==0)
		{
			BinOutput(loc+int2str(dumpnum)+".bin",sim,tess);
			++dumpnum;
		}
	  sim.TimeAdvance2Mid();
	  vector<int> removed=sim.RemoveCells(&remove);
	  sim.RefineCells(&refine,removed);
	  cntr++;
	}
      catch(UniversalError const& eo)
	{
	  cout << eo.GetErrorMessage() << endl;
	  cout << "Iteration number = " << sim.GetCycle() << endl;
	  for(int i=0;i<(int)eo.GetFields().size();++i)
	    {
	      cout << eo.GetFields()[i] << " = "
		   << eo.GetValues()[i] << endl;
	      if(eo.GetFields()[i]=="cell index")
		{
		  cout << "Exact value = " << eo.GetValues()[i] << endl;
		  cout << "Rounded values = " << (int)eo.GetValues()[i] << endl;
		  Vector2D temp(sim.GetMeshPoint((int)eo.GetValues()[i]));
		  cout << "cell x coordinate = " << temp.x << endl;
		  cout << "cell y coordinate = " << temp.y << endl;
		}
	    }
	  throw;
	}
      catch(vector<double> const& eo)
	{
	  cout << "Nan occurred in cell " << eo[0] <<  endl;
	  cout << "Cell density = " << eo[1] << endl;
	  cout << "Cell pressure = " << eo[2] << endl;
	  cout << "Cell x velocity = " << eo[3] << endl;
	  cout << "Cell y velocity = " << eo[4] << endl;
	  throw;
	}
      catch(char const* eo)
	{
	  cout << eo << endl;
	  throw;
	}
      catch(Primitive const& eo)
	{
	  for(int i=0;i<6;i++){
	    cout << eo[i] << endl;
	  }
	  throw;
	}
      // Infinite loop guard
      iter++;
      if(iter>max_iter)
	throw "Maximum number of iterations exceeded in main loop";
    }
  // Check if sim ran ok
  CheckSim(sim,tess);
  // Finalise
  return 0;
}


