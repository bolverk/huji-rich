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

  class HBCData
  {
  public:

    HBCData(Hllc const& rs,
	    IdealGas const& eos):
      hbc1_(rs),
      hbc3_(rs),
      hbc2_(CalcPrimitive(0.01,0.0015,
			  Vector2D(1,0),
			  eos),rs),
      hbc_(hbc2_,hbc1_,hbc3_,hbc3_) {}

    CustomOuter& getHBC(void)
    {
      return hbc_;
    }

  private:
    FreeFlow hbc1_;
    RigidWallHydro hbc3_;
    InFlow hbc2_;
    CustomOuter hbc_;
  };

  class RefineBow : public RefineStrategy
  {
  private:
    int N_;
  public:

    RefineBow(int npoints):
      N_(npoints) {}

    ~RefineBow(){}

    vector<int> CellsToRefine(Tessellation const& tess,
			      vector<Primitive> const& /*cells*/,
			      vector<vector<double> > const& /*tracers*/,
			      double /*time*/,
			      vector<Vector2D> &directions,
			      const vector<int> &Removed)
    {
      const double dv=1.2*5.0/(N_*N_);
      vector<int> res;
      directions.clear();
      Vector2D dir(1,0);
      for(int i=0;i<N_;++i)
	{
	  if(tess.GetVolume(i)>dv)
	    {
	      res.push_back(i);
	      directions.push_back(dir);
	    }
	}
      return RemoveDuplicatedLately(res,tess.GetPointNo(),directions,Removed);
    }
  };

  class RemoveBow : public RemovalStrategy
  {
  public:
    vector<int> CellsToRemove(Tessellation const& tess,
			      vector<Primitive> const& cells,
			      vector<vector<double> > const& /*tracers*/,
			      double /*time*/)const
    {
      vector<int> res;
      vector<double> merit;
      const int n=(int)cells.size();
      const double r=2/sqrt(1.0*n);
      for(int i=0;i<n;++i)
	{
	  const double dx=1-tess.GetMeshPoint(i).x;
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

  class SimData
  {
  public:

    SimData(double width=1, int np=40):
      tess_(),
      eos_(5./3.),
      outer_(-width, width,width*1.5,-width*1.5),
      rs_(),
      hbc_data_(rs_,eos_),
      interp_(eos_,outer_,hbc_data_.getHBC(), true, false, 0.1, 0.4),
      naive_(),
      point_motion_(naive_,hbc_data_.getHBC(),0.65,0.05,true,np,&outer_),
      force_(),
      refine_(np),
      remove_(),
      sim_(SquareMesh(np,np,width*2,width*3),
	   tess_,
	   interp_,
	   Uniform2D(0.01),
	   Uniform2D(0.0015),
	   Uniform2D(1),
	   Uniform2D(0),
	   eos_,
	   rs_,
	   point_motion_,
	   force_,
	   outer_,
	   hbc_data_.getHBC()) {}

    Tessellation& getTessellation(void)
    {
      return tess_;
    }

    RefineBow& getRefineScheme(void)
    {
      return refine_;
    }

    RemoveBow& getRemoveScheme(void)
    {
      return remove_;
    }

    hdsim& getSim(void)
    {
      return sim_;
    }

  private:
    VoronoiMesh tess_;
    const IdealGas eos_;
    const SquareBox outer_;
    const Hllc rs_;
    HBCData hbc_data_;
    LinearGaussConsistent interp_;
    Lagrangian naive_;
    RoundCells point_motion_;
    ZeroForce force_;
    RefineBow refine_;
    RemoveBow remove_;
    hdsim sim_;
  };

  void CheckSim(SimData& sim_data)
  {
    Tessellation const& tess = sim_data.getTessellation();
    hdsim const& sim = sim_data.getSim();
    const double GoodMass=0.06;
    const int n=tess.GetPointNo();
    double SimMass=0;
    for(int i=0;i<n;++i)
      SimMass+=tess.GetVolume(i)*sim.GetCell(i).Density;

    ofstream myfile("result.txt");
    myfile<<(GoodMass-SimMass)/GoodMass;
    myfile.close();
  }

  void main_loop(SimData& sim_data)
  {
    hdsim& sim = sim_data.getSim();

    sim.SetCfl(0.8);
    const int max_iter = 5e6;
    const double tf = 2;
    sim.SetEndTime(tf);
    while(tf>sim.GetTime()){
      try{
	sim.TimeAdvance2Mid();
	sim.RefineCells(&sim_data.getRefineScheme(),
			sim.RemoveCells(&sim_data.getRemoveScheme()));
      }
      catch(UniversalError const& eo){
	DisplayError(eo);
      }

      if(sim.GetCycle()>max_iter)
	throw UniversalError
	  ("Maximum number of iterations exceeded in main loop");
    }
  }
}

int main(void)
{
  SimData sim_data;

  main_loop(sim_data);

  CheckSim(sim_data);

  return 0;
}
