#include "PolyIntersect.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

typedef boost::mt19937_64 base_generator_type;

using std::vector;

//! \brief Overlap of two polygons
class PolygonOverlap
{
private:
  base_generator_type gen;
public:

  PolygonOverlap(void);

  double PolyArea(vector<Vector2D> const& polygon);
  
  double polygon_overlap_area(vector<Vector2D> const& ch1,
			      vector<Vector2D> const& ch2,double R0,double R1);
};
