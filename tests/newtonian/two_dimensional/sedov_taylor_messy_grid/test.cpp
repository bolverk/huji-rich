#ifdef RICH_MPI
#include <mpi.h>
#endif
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "source/tessellation/geometry.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/interpolations/pcm2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/step2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/test_2d/random_pert.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/mpi/mpi_macro.hpp"
#include "source/mpi/MeshPointsMPI.hpp"

using namespace std;
using namespace simulation2d;

namespace {

#ifdef RICH_MPI
  vector<Vector2D> process_positions(const SquareBox& boundary)
  {
    const Vector2D lower_left = boundary.getBoundary().first;
    const Vector2D upper_right = boundary.getBoundary().second;
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

  class SimData
  {
  public:

    SimData(void):
      width_(1),
      outer_(0,width_,width_,0),
#ifdef RICH_MPI
      proc_tess_(process_positions(outer_),outer_),
      init_points_(distribute_grid
		   (proc_tess_,
		    CartesianGridGenerator
		    (30,30, outer_.getBoundary().first,
		     outer_.getBoundary().second))),
#else
      init_points_(cartesian_mesh(30,30,outer_.getBoundary().first,
				  outer_.getBoundary().second)),
#endif
      tess_(),
      interp_method_(),
      density_(1),
      pressure_(0,0.06,0,0.06,1e4,0.01),
      xvelocity_(0),
      yvelocity_(0),
      eos_(5./3.),
      point_motion_(),
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
	   point_motion_,
	   force_,
	   outer_,
	   hbc_) {}

    hdsim& getSim(void)
    {
      return sim_;
    }

  private:
    const double width_;
    const SquareBox outer_;
#ifdef RICH_MPI
    VoronoiMesh proc_tess_;
#endif
    const vector<Vector2D> init_points_;
    VoronoiMesh tess_;
    PCM2D interp_method_;
    const Uniform2D density_;
    const Step2D pressure_;
    const Uniform2D xvelocity_;
    const Uniform2D yvelocity_;
    const IdealGas eos_;
    Lagrangian point_motion_;
    const Hllc rs_;
    const RigidWallHydro hbc_;
    ZeroForce force_;
    hdsim sim_;
  };

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.04,1e6);
    WriteTime diag("time.txt");
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
  write_snapshot_to_hdf5(sim, "process_"+int2str(get_mpi_rank())+"_final.h5");
#else
  write_snapshot_to_hdf5(sim, "final.h5");
#endif


#ifdef RICH_MPI
  MPI_Finalize();
#endif

  return 0;
}

