#ifdef RICH_MPI
#include "source/mpi/MeshPointsMPI.hpp"
#endif
#include <boost/foreach.hpp>
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/tessellation/shape_2d.hpp"
#include "source/newtonian/test_2d/piecewise.hpp"
#include "source/newtonian/two_dimensional/physical_geometry.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"
#include <numeric>

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

  vector<ComputationalCell> calc_init_cond
  (const Tessellation& tess)
  {
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i)
	{
      const Triangle triangle(Vector2D(0.5,0.6),
			      Vector2D(0.7,0.5),
			      Vector2D(0.4,0.4));
      const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
      res[i].density = 1;
      res[i].pressure = triangle(r) ? 2 : 1;
      res[i].velocity = triangle(r) ? Vector2D(1,-1) : Vector2D(0,0);
      res[i].tracers.push_back(triangle(r) ? 1 : 0);
    }
    return res;
  }

  #ifdef RICH_MPI

  vector<Vector2D> process_positions(const SquareBox& boundary)
  {
    const Vector2D lower_left = boundary.getBoundary().first;
    const Vector2D upper_right = boundary.getBoundary().second;
	int ws=0;
	MPI_Comm_size(MPI_COMM_WORLD,&ws);
    return RandSquare(ws,lower_left.x,upper_right.x,lower_left.y,upper_right.y);
  }

#endif
  
  class SimData
  {
  public:

    SimData(void):
      pg_(Vector2D(0,0),Vector2D(0,1)),
      outer_(Vector2D(0,0), Vector2D(1,1)),
#ifdef RICH_MPI
	proc_tess_(process_positions(outer_),outer_),
		init_points_(SquareMeshM(30,30,proc_tess_,outer_.getBoundary().first,outer_.getBoundary().second)),
		tess_(proc_tess_,init_points_,outer_),
#else
      init_points_(cartesian_mesh(30,30,
				  outer_.getBoundary().first,
				  outer_.getBoundary().second)),
		tess_(init_points_,outer_),
#endif

      triangle_(Vector2D(0.5,0.6),
		Vector2D(0.7,0.5),
		Vector2D(0.4,0.4)),
      eos_(5./3.),
      pm_(),
      evc_(),
      rs_(),
      force_(pg_.getAxis()),
      tsf_(0.3),
      fc_(rs_),
      eu_(),
      cu_(),
      sim_(
	  #ifdef RICH_MPI
		  proc_tess_,
#endif
	tess_,
	   outer_,
	   pg_,
	   calc_init_cond(tess_),
	   eos_,
	   pm_,
	   evc_,
	   force_,
	   tsf_,
	   fc_,
	   eu_,
	   cu_,
	   TracerStickerNames (vector<string> (1,"tracer"),vector<string>())) {}

    hdsim& getSim(void)
    {
      return sim_;
    }

  private:
    CylindricalSymmetry pg_;
    SquareBox outer_;
#ifdef RICH_MPI
    VoronoiMesh proc_tess_;
#endif
    const vector<Vector2D> init_points_;
    VoronoiMesh tess_;
    const Triangle triangle_;
    IdealGas eos_;
    Lagrangian pm_;
    const StationaryBox evc_;
    Hllc rs_;
    CylindricalComplementary force_;
    SimpleCFL tsf_;
    SimpleFluxCalculator fc_;
    const SimpleExtensiveUpdater eu_;
    SimpleCellUpdater cu_;
    hdsim sim_;
  };

  class WriteConserved: public DiagnosticFunction
  {
  public:

    WriteConserved(string const& fname):
      cons_(), fname_(fname)
#ifdef RICH_MPI
	,rank_(0)
#endif
	  {
#ifdef RICH_MPI
		MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
#endif
	  }

    void operator()(hdsim const& sim)
    {
      cons_.push_back(accumulate(next(sim.getAllExtensives().begin()),
				 sim.getAllExtensives().end(),
				 sim.getAllExtensives().front()));
    }

    ~WriteConserved(void)
    {
#ifdef RICH_MPI
      if(rank_==0){
#endif
	ofstream f(fname_.c_str());
	for(size_t i=0;i<cons_.size();++i)
	  f << cons_[i].mass << " "
	    << cons_[i].momentum.x << " "
	    << cons_[i].momentum.y << " "
	    << cons_[i].energy << " "
	    << cons_[i].tracers[0] << "\n";
	f.close();
#ifdef RICH_MPI
      }
#endif
    }

  private:
    mutable vector<Extensive> cons_;
    const string fname_;
#ifdef RICH_MPI
	int rank_;
#endif
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

