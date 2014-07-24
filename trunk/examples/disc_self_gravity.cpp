#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/Circle2D.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/source_terms/SelfGravity2D.hpp"
#include "source/newtonian/two_dimensional/custom_evolutions/ConstantPrimitiveEvolution.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/misc/int2str.hpp"
#include "source/tessellation/RoundGrid.hpp"

class KT_Density :public SpatialDistribution
{
private:
  double M_,ROuter_;
public:
  KT_Density(double M,double ROuter):M_(M),ROuter_(ROuter){};
  ~KT_Density(){};
  double operator()(Vector2D const& point) const
  {
    if(abs(point)<ROuter_)
      return M_*pow(1+abs(point)*abs(point),-1.5)/(2*3.1415);
    else
      return 0.01*M_*pow(1+abs(point)*abs(point),-1.5)/(2*3.1415);
  }
};

class KT_xvel :public SpatialDistribution
{
private:
  double M_;
public:
  KT_xvel(double M):M_(M){};
  ~KT_xvel(){};
  double operator()(Vector2D const& point) const
  {
    return -sqrt(M_*pow(1+abs(point)*abs(point),-1.5))*point.y;
  }
};

class KT_yvel :public SpatialDistribution
{
private:
  double M_;
public:
  KT_yvel(double M):M_(M){};
  ~KT_yvel(){};
  double operator()(Vector2D const& point) const
  {
    return sqrt(M_*pow(1+abs(point)*abs(point),-1.5))*point.x;
  }
};

int main(void)
{
  // Set up the initial grid points
  int np = 50;
  double InnerR=0.1;
  double OuterR=10;
  vector<Vector2D> InitPoints=CirclePointsRmax_a
    (2*np*np,
     InnerR*1.35,OuterR,0,
     0,OuterR,OuterR,-OuterR,-OuterR,-0.6);
  vector<Vector2D> OuterCircle1 = 
    circle_circumference(np*2,
			 OuterR*(1+0.5/np),
			 Vector2D(0,0));
  vector<Vector2D> OuterCircle2 = 
    circle_circumference(np*2,
			 OuterR*(1+0.25/np),
			 Vector2D(0,0));
  vector<Vector2D> OuterCircle3 = 
    circle_circumference(np*2,
			 OuterR*(1+0.75/np),
			 Vector2D(0,0));
  /*
    vector<Vector2D> OuterCircle1=Circle(np*2,OuterR*(1+0.5/np));
    vector<Vector2D> OuterCircle2=Circle(np*2,OuterR*(1+1.25/np));
    vector<Vector2D> OuterCircle3=Circle(np*2,OuterR*(1+1.75/np));
  */
  vector<Vector2D> InnerCircle=CirclePointsRmax(np,0.1*InnerR,InnerR);
  InitPoints.insert(InitPoints.begin(),InnerCircle.begin(),InnerCircle.end());
  InitPoints.insert(InitPoints.begin(),OuterCircle1.begin(),OuterCircle1.end());
  InitPoints.insert(InitPoints.begin(),OuterCircle2.begin(),OuterCircle2.end());
  InitPoints.insert(InitPoints.begin(),OuterCircle3.begin(),OuterCircle3.end());

  // Set up the boundary type for the points
  double BoxOuter=OuterR*(1+2.25/np);
  SquareBox outer(-BoxOuter,BoxOuter,BoxOuter,-BoxOuter);

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
  RoundCells pointmotion(l_motion,hbc,1,0.1,6*np);

  // Set up the interpolation
  LinearGaussConsistent interpolation(eos,outer,hbc);

  // Set up the initial Hydro
  double M=1;
  double p=0.02;
  Circle2D pressure(0,0,OuterR*(1+1.5/np),p,p*0.001);
  KT_Density density(M,OuterR*(1+1.5/np));
  KT_xvel xvelocity(M);
  KT_yvel yvelocity(M);

  // Set up the external source term
  SelfGravity acc(0.45,0.1);
  ConservativeForce force(acc);

  // Set up the simulation
  hdsim sim(InitPoints,tess,interpolation,density,pressure,xvelocity,
	    yvelocity,eos,rs,pointmotion,force,outer,hbc,false,false);
	
  // Give custom evolution to edge
  ConstantPrimitiveEvolution cp;
  sim.custom_evolution_manager.addCustomEvolution(&cp);
  for(int i=0;i<4*np;++i)
    sim.custom_evolution_indices[i] = 1;

  // Choose the Courant number
  double cfl=0.7;
  sim.SetCfl(cfl);

  // How long shall we run the simulation?
  double tend=25;
  sim.SetEndTime(tend);

  // Custom output criteria
  double output_dt=0.3;
  double last_dump_time=0;
  int dump_number=0;

  write_snapshot_to_hdf5(sim,"c:\\sim_data\\output0.bin");

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
