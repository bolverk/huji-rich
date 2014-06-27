#ifdef RICH_MPI
#include "../../mpi/mpi_macro.hpp"
#include "../../mpi/MeshPointsMPI.hpp"
#endif
#include "triangle_step.hpp"
#include "../../misc/universal_error.hpp"

namespace {

  #ifdef RICH_MPI
  class PointsOnCircle: public Index2Member<Vector2D>
  {
  public:

    PointsOnCircle(const Vector2D& center,
		   double radius,
		   int n):
      center_(center), radius_(radius), n_(n) {}

    size_t getLength(void) const
    {
      return (size_t)n_;
    }

    Vector2D operator()(size_t i) const
    {
      if(n_>3){
	const double angle = 2*M_PI*(double)(i-1)/(double)(n_-1);
	return radius_*Vector2D(cos(angle), sin(angle)) + center_;
      }
      else{
	if(i==0)
	  return center_;
	else{
	  const double angle = 2*M_PI*(double)(i-1)/(double)(n_-1);
	  return radius_*Vector2D(cos(angle), sin(angle)) + center_;
	}
      }
    }

  private:
    const Vector2D& center_;
    const double radius_;
    const int n_;
  };

  vector<Vector2D> process_positions(const SquareBox& boundary)
  {
    const Vector2D lower_left = boundary.getBoundary().first;
    const Vector2D upper_right = boundary.getBoundary().second;
    /*
    const Vector2D center = 0.5*(lower_left+upper_right);
    const double radius = 0.25*min((upper_right-lower_left).x,
				   (upper_right-lower_left).y);
    return serial_generate(PointsOnCircle(center,
					  radius,
					  get_mpi_size()));
    */
    vector<Vector2D> res(get_mpi_size());
    if(get_mpi_rank()==0){
      res = RandSquare(get_mpi_size(),
		       lower_left.x,upper_right.x,
		       lower_left.y,upper_right.y);
    }
    MPI_VectorBcast_Vector2D(res,0,MPI_COMM_WORLD,get_mpi_rank());
    return res;
  }
  #endif

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

TriangleStep::TriangleStep(string const& point_motion, double width):
  outer_(0,width,width,0),
  #ifdef RICH_MPI
  proc_tess_(process_positions(outer_),outer_),
  init_points_(distribute_grid
	       (proc_tess_,
		CartesianGridGenerator(30,30,
				       outer_.getBoundary().first,
				       outer_.getBoundary().second))),
  #else
  init_points_(cartesian_mesh(30,30,
			      outer_.getBoundary().first,
			      outer_.getBoundary().second)),
  #endif
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
       #ifdef RICH_MPI
       proc_tess_,
       #endif
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

TriangleStep::~TriangleStep(void) {}
