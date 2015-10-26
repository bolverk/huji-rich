#include <iostream>
#include "source/tessellation/geometry.hpp"
#include "source/tessellation/Delaunay.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/misc/mesh_generator.hpp"
#ifdef RICH_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#endif // RICH_MPI
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include <boost/foreach.hpp>
#include "source/tessellation/VoronoiMesh.hpp"

using namespace std;

#ifdef RICH_MPI
vector<Vector2D> delineate_rectangle
(const pair<Vector2D,Vector2D>& boundary)
{
  vector<Vector2D> res(4);
  res.at(0) = boundary.first;
  res.at(1) = Vector2D(boundary.second.x,boundary.first.y);
  res.at(2) = boundary.second;
  res.at(3) = Vector2D(boundary.first.x,boundary.second.y);
  return res;
}

vector<Vector2D> my_convex_hull
(const Tessellation& tess,
 int index)
{
  vector<Vector2D> res;
  ConvexHull(res,tess,index);
  return res;
}

vector<Vector2D> distribute_grid
(const vector<Vector2D>& complete_grid,
 const Tessellation& proc_tess)
{
  const boost::mpi::communicator world;
  vector<Vector2D> res;
  const vector<Vector2D> ch_list =
    my_convex_hull(proc_tess,world.rank());
  BOOST_FOREACH(const Vector2D& v, complete_grid)
    {
      if(PointInCell(ch_list,v))
	res.push_back(v);
    }
  return res;
}
#endif // RICH_MPI

int main(void)
{
#ifdef RICH_MPI

  boost::mpi::environment env;
  boost::mpi::communicator world;
  
  try{

  const SquareBox obc
    (Vector2D(0,0),
     Vector2D(1,1));
  const vector<Vector2D> all_points = 
    RandSquare(1000,
	       obc.getBoundary().first.x,
	       obc.getBoundary().second.x,
	       obc.getBoundary().first.y,
	       obc.getBoundary().second.y);
  if(world.rank()==0){
    Delaunay tri;
    tri.build_delaunay
      (all_points, 
       delineate_rectangle
       (obc.getBoundary()));
    tri.BuildBoundary(&obc, obc.GetBoxEdges());
    WriteDelaunay(tri,string("whole.h5"));
  }
  VoronoiMesh super
    (RandSquare(world.size(),
		obc.getBoundary().first.x,
		obc.getBoundary().second.x,
		obc.getBoundary().first.y,
		obc.getBoundary().second.y),
     obc);
  const vector<Vector2D> local_points =
    distribute_grid
    (all_points,
     super);
  Delaunay tri;
  tri.build_delaunay
    (local_points,
     delineate_rectangle
     (obc.getBoundary()));
  vector<vector<int> > dummy_1;
  //vector<int> dummy_2;
  tri.BuildBoundary
    (&obc,
     super,
     dummy_1);
  WriteDelaunay(tri,string("part_"+int2str(world.rank())+".h5"));
  }
  catch(const UniversalError& eo){
    cout << eo.GetErrorMessage() << endl;
    throw;
  }
#else

  write_number(0,"serial_ignore.txt");

#endif // RICH_MPI
  return 0;
}
