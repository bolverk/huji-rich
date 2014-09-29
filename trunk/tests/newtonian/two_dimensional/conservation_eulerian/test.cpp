#ifdef RICH_MPI
#include <mpi.h>
#endif
#include <boost/foreach.hpp>
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/tessellation/shape_2d.hpp"
#include "source/newtonian/test_2d/piecewise.hpp"
#include "source/newtonian/two_dimensional/physical_geometry.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/interpolations/pcm2d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"

using namespace std;
using namespace simulation2d;

namespace {

  bool same_side(const Vector2D& p1,
		 const Vector2D& p2,
		 const Vector2D& a,
		 const Vector2D& b)
  {
    const double cp1 = CrossProduct(b-a,p1-a);
    const double cp2 = CrossProduct(b-a,p2-a);
    return cp1*cp2>=0;
  }

  class Triangle: public Shape2D
  {
  public:

    const Vector2D p1;
    const Vector2D p2;
    const Vector2D p3;

    Triangle(const Vector2D& p1_i,
	     const Vector2D& p2_i,
	     const Vector2D& p3_i):
      p1(p1_i), p2(p2_i), p3(p3_i) {}

    bool operator()(const Vector2D& p) const
    {
      return same_side(p,p1,p2,p3)&&
	same_side(p,p2,p1,p3)&&
	same_side(p,p3,p1,p2);
    }
  };

  class SimData
  {
  public:

    SimData(void):
      outer_(Vector2D(0,0), Vector2D(1,1)),
#ifdef RICH_MPI
      proc_tess_(process_positions(outer_),outer_),
      init_points_
      (distribute_grid(proc_tess_,CartesianGridGenerator
		       (30,30, outer_.getBoundary().first,
			outer_.getBoundary().second))),
#else
      init_points_(cartesian_mesh(30,30,
				  outer_.getBoundary().first,
				  outer_.getBoundary().second)),
#endif
      tess_(),
      interp_method_(),
      triangle_(Vector2D(0.5,0.6),
		Vector2D(0.7,0.5),
		Vector2D(0.4,0.4)),
      eos_(5./3.),
      pm_(),
      rs_(),
      hbc_(rs_),
      force_(),
      sim_(init_points_,
	   tess_,
	   #ifdef RICH_MPI
	   proc_tess_,
	   #endif
	   interp_method_,
	   Uniform2D(1),
	   Piecewise(triangle_,
		     Uniform2D(2),
		     Uniform2D(1)),
	   Piecewise(triangle_,
		     Uniform2D(1),
		     Uniform2D(0)),
	   Piecewise(triangle_,
		     Uniform2D(-1),
		     Uniform2D(0)),
	   eos_,
	   rs_,
	   pm_,
	   force_,
	   outer_,
	   hbc_) {}

    hdsim& getSim(void)
    {
      return sim_;
    }

  private:
    SquareBox outer_;
#ifdef RICH_MPI
    VoronoiMesh proc_tess_;
#endif
    const vector<Vector2D> init_points_;
    VoronoiMesh tess_;
    PCM2D interp_method_;
    const Triangle triangle_;
    IdealGas eos_;
    Eulerian pm_;
    Hllc rs_;
    RigidWallHydro hbc_;
    ZeroForce force_;
    hdsim sim_;
  };

  class WriteConserved: public DiagnosticFunction
  {
  public:

    WriteConserved(string const& fname):
      cons_(), fname_(fname) {}

    void operator()(hdsim const& sim)
    {
      cons_.push_back(total_conserved(sim));
    }

    ~WriteConserved(void)
    {
#ifdef RICH_MPI
      if(get_mpi_rank()==0){
#endif
	ofstream f(fname_.c_str());
	for(size_t i=0;i<cons_.size();++i)
	  f << cons_[i].Mass << " "
	    << cons_[i].Momentum.x << " "
	    << cons_[i].Momentum.y << " "
	    << cons_[i].Energy << "\n";
	f.close();
#ifdef RICH_MPI
      }
#endif
    }

  private:
    mutable vector<Conserved> cons_;
    const string fname_;
  };
}

namespace {
  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.05,1e6);
    WriteConserved diag("res.txt");
    main_loop(sim, 
	      term_cond, 
	      &hdsim::TimeAdvance,
	      &diag);
  }
}

int main(void)
{
#ifdef RICH_MPI
  MPI_Init(NULL, NULL);
#endif
  SimData sim_data;
  hdsim& sim = sim_data.getSim();

  my_main_loop(sim);

#ifdef RICH_MPI
  MPI_Finalize();
#endif

  return 0;
}

