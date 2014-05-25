#include <cstdlib>
#include "random_pert.hpp"

vector<Vector2D> add_random_pert(vector<Vector2D> const& original,
				 double amp_x, double amp_y)
{
  srand(0);
  vector<Vector2D> res = original;
  for(int i=0;i<(int)res.size();++i){
    res[i].x += 
      amp_x*2*((double)rand()/(double)RAND_MAX-0.5);
    res[i].y += 
      amp_y*2*((double)rand()/(double)RAND_MAX-0.5);
  }
  return res;
}
