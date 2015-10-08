#include <iostream>
#include "source/tessellation/geometry.hpp"
#include "source/tessellation/Delaunay.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/misc/mesh_generator.hpp"

using namespace std;

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

int main(void)
{
#ifdef RICH_MPI
  const SquareBox obc
    (Vector2D(0,0),
     Vector2D(1,1));
  const vector<Vector2D> all_points = 
    RandSquare(1000,
	       obc.getBoundary().first.x,
	       obc.getBoundary().second.x,
	       obc.getBoundary().first.y,
	       obc.getBoundary().second.y);
    Delaunay tri;
    tri.build_delaunay
      (all_points, 
       delineate_rectangle
       (obc.getBoundary()));
    
#else

  write_number(0,"serial_ignore.txt");

#endif // RICH_MPI
  return 0;
}
