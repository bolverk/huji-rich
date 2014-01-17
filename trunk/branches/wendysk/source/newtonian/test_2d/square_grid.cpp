#include <cmath>
#include "square_grid.hpp"

vector<Vector2D> square_grid(double side, int np)
{
  vector<Vector2D> res((int)pow((double)np,2));
  const double width = side/(double)np;
  for(int i=0;i<np;++i){
    for(int j=0;j<np;++j){
      Vector2D point;
      point.x = ((double)i+0.5)*width;
      point.y = ((double)j+0.5)*width;
      res[i*np+j] = point;
    }
  }
  return res;
}

vector<Vector2D> offset_grid(vector<Vector2D> const& original,
			     Vector2D const& offset)
{
  vector<Vector2D> res(original.size());
  for(int i=0;i<(int)original.size();++i)
    res[i] = original[i] + offset;
  return res;
}
