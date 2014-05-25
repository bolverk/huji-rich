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
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/interpolations/linear_gauss.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/InnerRigid.hpp"
#include "source/newtonian/two_dimensional/custom_evolutions/RigidBodyEvolve.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/FreeFlow.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/CustomOuter.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/InFlow.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"

using namespace std;

namespace {
vector<Vector2D> CustomPoints(int nx,int ny,double sidex,double sidey)
{
  vector<Vector2D> res;
  double widthx = sidex/(double)nx;
  double widthy = sidey/(double)ny;
  Vector2D point;
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      point.x = ((double)i+0.5)*widthx-sidex/2;
      point.y = ((double)j+0.5)*widthy-sidey/2;
      if(point.y>0.1*point.x-0.428)
	res.push_back(point);
    }
  }
  return res;
}

vector<Vector2D> Line(int PointNum,double xmin,double xmax,double ymin,double ymax)
{
  double dy=ymax-ymin,dx=xmax-xmin;
  double angle=atan2(dy,dx);
  double length=sqrt(dy*dy+dx*dx);
  double dl=length/PointNum;

  vector<Vector2D> res;
  Vector2D temp;

  for(int i=0;i<PointNum;++i)
    {
      temp.Set(xmin+dl*(i+0.5)*cos(angle),ymin+dl*(i+0.5)*sin(angle));
      res.push_back(temp);
    }
  return res;
}

void write_output(hdsim const& sim)
{
  write_generating_points(sim,"mesh_points.txt");
  write_cells_property(sim,"velocity y","yvelocities.txt");
}
}

int main(void)
{
  // Initialization
  double width = 1;
  int np = 30;
  int innerNum=np*2;

  // Define Outer Box
  SquareBox outer(-width/2,width/2,width/2,-width/2);

  // Create mesh points
  vector<Vector2D> InitPoints=CustomPoints(np,np,width,width);
  vector<Vector2D> InnerBoundary1=Line(innerNum,-width/2+0.01,width/2,-width/2,width/2-0.9);
  vector<Vector2D> InnerBoundary2;
  Vector2D temp;
  double alpha=atan((InnerBoundary1[0].y-InnerBoundary1[1].y)/
		    (InnerBoundary1[0].x-InnerBoundary1[1].x));
  write_number(alpha,"wedge_angle.txt");
  double a=0.01;
  for(int i=0;i<(int)InnerBoundary1.size();++i)
    {
      temp.Set(InnerBoundary1[i].x-a*sin(alpha),
	       InnerBoundary1[i].y+a*cos(alpha));
      InnerBoundary2.push_back(temp);
    }
	
  InitPoints.insert(InitPoints.begin(),InnerBoundary2.begin(),InnerBoundary2.end());
  InitPoints.insert(InitPoints.begin(),InnerBoundary1.begin(),InnerBoundary1.end());

  //EOS
  IdealGas eos(5./3.);
  write_number(eos.getAdiabaticIndex(),"adiabatic_index.txt");
  Hllc rs;

  // Hydro boundary conditions
  FreeFlow f_hbc(rs);
  RigidWallHydro r_hbc(rs);

  Primitive InFlux;
  InFlux.Density=1;
  InFlux.Pressure=1;
  InFlux.Velocity.Set(4,0);
  InFlux.Energy=eos.dp2e(InFlux.Density,InFlux.Pressure);
  InFlux.SoundSpeed=eos.dp2c(InFlux.Density,InFlux.Pressure);
  write_number(InFlux.Velocity.x/InFlux.SoundSpeed,
	       "mach_number.txt");

  InFlow i_hbc(InFlux, rs);

  CustomOuter hbc(&i_hbc,&f_hbc,&r_hbc,&f_hbc);	

  // Create tessellation structure
  VoronoiMesh tess;

  // Spatial interpolation
  LinearGauss pcm2d(outer,&hbc,true,false);

  // Hydro initial conditions

  Uniform2D pressure(1);
  Uniform2D density(1);
  Uniform2D xvelocity(4);
  Uniform2D yvelocity(0);
	
  //Point Motion
  Eulerian pointmotion;

  //External Forces
  ZeroForce force;

  hdsim sim(InitPoints, //Initial points
	    &tess, // Tessellation
	    &pcm2d, // Spatial reconstruction
	    density, pressure, xvelocity, yvelocity,
	    eos,
	    rs,
	    &pointmotion,
	    &force,
	    &outer,
	    &hbc);

  // Assign the inner boundary
  RigidBodyEvolve rigid;
  for(int i=0;i<innerNum;++i)
    sim.CellsEvolve[i]=&rigid;

  // Main process
  sim.SetCfl(0.5);
  //  int iter = 0;
  const int max_iter = (int) 5e6;

	
  double tf =1;
  sim.SetEndTime(tf);
	
  //Main simulation run
  while(tf>sim.GetTime()){
    try{
      sim.TimeAdvance2Mid();
    }
    catch(UniversalError const& eo){
      cout << eo.GetErrorMessage() << endl;
      for(int i=0;i<(int)eo.GetFields().size();++i){
	cout << eo.GetFields()[i] << " = "
	     << eo.GetValues()[i] << endl;
	cout << "Iteration number = " << sim.GetCycle() << endl;
	if(eo.GetFields()[i]=="cell index"){
	  cout << "Exact value = " << eo.GetValues()[i] << endl;
	  cout << "Rounded values = " << (int)eo.GetValues()[i] << endl;
	  Vector2D temp2(sim.GetMeshPoint((int)eo.GetValues()[i]));
	  cout << "cell x coordinate = " << temp2.x << endl;
	  cout << "cell y coordinate = " << temp2.y << endl;
	}
      }
      throw;
    }
    catch(vector<double> const& eo){
      cout << "Nan occurred in cell " << eo[0] <<  endl;
      cout << "Cell density = " << eo[1] << endl;
      cout << "Cell pressure = " << eo[2] << endl;
      cout << "Cell x velocity = " << eo[3] << endl;
      cout << "Cell y velocity = " << eo[4] << endl;
      throw;
    }
    catch(char const* eo){
      cout << eo << endl;
      throw;
    }
    catch(Primitive const& eo){
      for(int i=0;i<6;i++){
	cout << eo[i] << endl;
      }
      throw;
    }

    // Infinite loop guard
    if(sim.GetCycle()>max_iter)
      throw "Maximum number of iterations exceeded in main loop";
  }

  write_output(sim);

  return 0;
}


