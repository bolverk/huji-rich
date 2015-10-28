#ifdef RICH_MPI
#include "source/mpi/MeshPointsMPI.hpp"
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
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/tessellation/right_rectangle.hpp"
#include "source/newtonian/test_2d/piecewise.hpp"
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include <fenv.h>
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"

using namespace std;
using namespace simulation2d;

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

namespace {
  vector<ComputationalCell> calc_init_cond(const Tessellation& tess)
  {
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i){
      res[i].density = 1;
      res[i].pressure = tess.GetMeshPoint(static_cast<int>(i)).x < 0.5 ?
								   2 : 1;
      res[i].velocity = Vector2D(0,0);
    }
    return res;
  }
}

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
		  (30,30,outer_.getBoundary().first,
		   outer_.getBoundary().second))),
    #else
    pg_(),
    init_points_(cartesian_mesh(30,30,outer_.getBoundary().first,
				outer_.getBoundary().second)),
    #endif
    tess_(init_points_,outer_),
    eos_(5./3.),
    pm_naive_(),
    rs_(),
    point_motion_(pm_naive_,eos_),
    evc_(),
    force_(),
    tsf_(0.3),
    fc_(rs_),
    eu_(),
    cu_(),
    sim_(tess_,
	 #ifdef RICH_MPI
	 proc_tess_,
	 #endif
	 outer_,
	 pg_,
	 calc_init_cond(tess_),
	 eos_,
	 point_motion_,
	 evc_,
	 force_,
	 tsf_,
	 fc_,
	 eu_,
	 cu_) {}

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
  const SlabSymmetry pg_;
  const vector<Vector2D> init_points_;
  VoronoiMesh tess_;
  const IdealGas eos_;
  Lagrangian pm_naive_;
  const Hllc rs_;
  RoundCells point_motion_;
  const StationaryBox evc_;
  ZeroForce force_;
  const SimpleCFL tsf_;
  const SimpleFluxCalculator fc_;
  const SimpleExtensiveUpdater eu_;
  const SimpleCellUpdater cu_;
  hdsim sim_;
};

namespace {
  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.225, 1e6);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      &hdsim::TimeAdvance,
	      &diag);
  }
}

int main(void)
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  #ifdef RICH_MPI
  MPI_Init(NULL, NULL);
  #endif

  SimData sim_data;
  hdsim& sim = sim_data.getSim();

  my_main_loop(sim);

  #ifdef RICH_MPI
  write_snapshot_to_hdf5(sim,"process_"+int2str(get_mpi_rank())+"_final.h5");
  #else
  write_snapshot_to_hdf5(sim,"final.h5");
  #endif

  #ifdef RICH_MPI
  MPI_Finalize();
  #endif

  return 0;
}

