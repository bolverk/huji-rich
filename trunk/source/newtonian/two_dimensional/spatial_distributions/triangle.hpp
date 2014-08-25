/*! \file triangle.hpp
  \brief Triangle step function
  \author Almog Yalinewich
*/

#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP 1

#include "../spatial_distribution2d.hpp"
#include "../../../tessellation/geometry.hpp"

//! \brief A spatial distribution with a triangle step function
class Triangle: public SpatialDistribution
{
public:

  /*! \brief Class constructor
    \param vv Vector of vertices
    \param vi Value inside the triangle
    \param vo Value outside the triangle
   */
  Triangle(vector<Vector2D> vv, double vi, double vo);

  double operator()(Vector2D const& r) const;

private:

  vector<Vector2D> vv_;
  double vi_;
  double vo_;
};

#endif // TRIANGLE_HPP