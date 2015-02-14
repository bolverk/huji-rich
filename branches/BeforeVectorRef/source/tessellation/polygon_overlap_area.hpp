/*! \file polygon_overlap_area.hpp
  \brief Calculates the overlap between two polygons
  \author Elad Steinberg
 */

#include "PolyIntersect.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

//! \brief Alias for random number generator
typedef boost::mt19937_64 base_generator_type;

using std::vector;

//! \brief Overlap of two polygons
class PolygonOverlap
{
private:
  base_generator_type gen;
public:
//! \brief Class constructor
  PolygonOverlap(void);

  /*!
  \brief Calcualted the area of a convex polygon
  \param polygon The vertices of the polygon, should be ordered as a convex hull
  \return The area
  */
  double PolyArea(vector<Vector2D> const& polygon);
  /*!
  \brief Calculates the area overlaped between two convex polygons
  \param ch1 The convex hull of the first polygon
  \param ch2 The convex hull of the second polygon
  \param R0 A fraction of the effective radius of the first polygon 
  \param R1 A fraction of the effective radius of the second polygon
  \return The overlaped area
  */
  double polygon_overlap_area(vector<Vector2D> const& ch1,
			      vector<Vector2D> const& ch2,double R0,double R1);
};
