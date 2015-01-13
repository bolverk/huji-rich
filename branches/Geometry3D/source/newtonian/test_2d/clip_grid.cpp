#include <boost/foreach.hpp>
#include "clip_grid.hpp"

vector<Vector2D> clip_grid(const Shape2D& shape,
			   const vector<Vector2D>& original)
{
  vector<Vector2D> res;
  BOOST_FOREACH(Vector2D point, original)
    {
      if(shape(point))
	res.push_back(point);
    }
  return res;
}
