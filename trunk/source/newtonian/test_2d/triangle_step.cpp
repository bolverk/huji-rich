#include "triangle_step.hpp"
#include "../../misc/universal_error.hpp"

namespace {
vector<Vector2D> get_triangle_vertices(void)
{
  vector<Vector2D> res;
  res.push_back(Vector2D(0.5,0.6));
  res.push_back(Vector2D(0.7,0.5));
  res.push_back(Vector2D(0.4,0.4));
  return res;
}

PointMotion& choose_between
(string const& choice,
 string const& name1, PointMotion& opt1,
 string const& name2, PointMotion& opt2)
{
  if(name1==choice)
    return opt1;
  else if(name2==choice)
    return opt2;
  else
    throw UniversalError("Unknown option "+choice);
}
}

TriangleStep::TriangleStep(string const& point_motion):
  width_(1),
  init_points_(square_grid(width_,30)),
  outer_(0,width_,width_,0),
  tess_(),
  interp_method_(),
  density_(1),
  triangle_vertices_(get_triangle_vertices()),
  pressure_(triangle_vertices_,2,1),
  xvelocity_(triangle_vertices_,1,0),
  yvelocity_(triangle_vertices_,-1,0),
  eos_(5./3.),
  eulerian_(),
  lagrangian_(),
  rs_(),
  hbc_(rs_),
  force_(),
  sim_(init_points_,
       tess_,
       interp_method_,
       density_,
       pressure_,
       xvelocity_,
       yvelocity_,
       eos_,
       rs_,
       choose_between(point_motion,
		      "eulerian",
		      eulerian_,
		      "lagrangian",
		      lagrangian_),
       force_,
       outer_,
       hbc_) {}

hdsim& TriangleStep::getSim(void)
{
  return sim_;
}
