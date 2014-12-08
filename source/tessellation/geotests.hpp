/*! \file geotests.hpp
  \brief Various checks for geometric data
  \author Elad Steinberg
 */
#ifndef GEOTEST_HPP
#define GEOTEST_HPP 1

#include <vector>
#include "geometry.hpp"
#include "exactmath.hpp"

/*! \brief Checks whether 3 given points are on a counterclockwise circle, clockwise circle or colinear.
 \param points The three points to check
 \return The result will be positive, negative or 0 respectively.
 */
double orient2d(boost::array<Vector2D,3> const& points);
/*! \brief Checks whether 3 given points are on a counterclockwise circle, clockwise circle or colinear using adpative math.
 \param points The three points to check
 \param detsum The determinant obtained from orient2d
 \return The result will be positive, negative or 0 respectively.
 */
double orient2dAdapt(boost::array<Vector2D,3> const& points, double detsum);

/*!\brief Checks whether the 4th point is inside, outside or on the counterclockwise circle created by the first 3 points.
  \param points The points to check
  \return The result will be positive, negative or 0 respectively.
 */
double incircle(const Vector2D& point_1,
		const Vector2D& point_2,
		const Vector2D& point_3,
		const Vector2D& point_4);

/*!\brief Checks whether the 4th point is inside, outside or on the counterclockwise circle created by the first 3 points using adaptive math.
  \param points The points to check
  \param permanent The error estimate
  \return The result will be positive, negative or 0 respectively.
 */
double incircleadapt(const Vector2D& point_1,
		     const Vector2D& point_2,
		     const Vector2D& point_3,
		     const Vector2D& point_4,
		     double permanent);

#endif // GEOTEST_HPP
