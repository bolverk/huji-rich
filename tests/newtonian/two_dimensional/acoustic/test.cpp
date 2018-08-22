#include <iostream>
#include "source/newtonian/test_1d/acoustic.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/PeriodicBox.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/test_2d/profile_1d.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/modular_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/ghost_point_generators/PeriodicGhostGenerator.hpp"
#include "source/newtonian/two_dimensional/interpolations/LinearGaussImproved.hpp"
#include "source/newtonian/two_dimensional/periodic_edge_velocities.hpp"
#include "source/tessellation/RoundGrid.hpp"
#ifdef RICH_MPI
#include "source/mpi/MeshPointsMPI.hpp"
#endif // RICH_MPI

using namespace std;
using namespace simulation2d;

namespace {

  void report_error(UniversalError const& eo)
  {
    cout << "Caught universal error" << endl;
    cout << eo.GetErrorMessage() << endl;
    for(size_t i = 0;i<eo.GetFields().size();++i){
      cout << eo.GetFields()[i] << " = "
	   << eo.GetValues()[i] << endl;
    }
  }

  vector<ComputationalCell> calc_init_cond(const Tessellation& tess,
					   const EquationOfState& eos,
					   double width)
  {
    const AcousticInitCond aic
      (read_number("ambient_density.txt"),
       read_number("ambient_pressure.txt"),
       eos,
       read_number("amplitude.txt"),
       width);
    const SpatialDistribution1D& density_dist = aic.getProfile("density");
    const SpatialDistribution1D& pressure_dist = aic.getProfile("pressure");
    const SpatialDistribution1D& velocity_dist = aic.getProfile("xvelocity");
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i){
      const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
      res[i].density = density_dist(r.x);
      res[i].pressure = pressure_dist(r.x);
      res[i].velocity = Vector2D(velocity_dist(r.x),0);
    }
    return res;
  }

#ifdef RICH_MPI

  vector<Vector2D> create_meta_points(const PeriodicBox& outer)
  {
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    return RoundGrid(RandSquare(world_size,
				outer.getBoundaries().first,
				outer.getBoundaries().second),
		     &outer);
  }

  vector<Vector2D> create_mesh_generating_points(int np,
						 const  PeriodicBox& outer,
						 const Tessellation& meta)
  {
    return RoundGrid(RandSquare(2*np*np,
				meta,
				outer.getBoundaries().first,
				outer.getBoundaries().second),
		     &outer, 30, &meta);
  }

#else
  /*

  vector<Vector2D> create_mesh_generating_points
  (int np,
   const PeriodicBox& outer)
  {
    vector<Vector2D> res = RandSquare
      (2*np*np,
       outer.getBoundaries().first,
       outer.getBoundaries().second);
    res = RoundGrid(res,&outer,30);
    return res;
  }
  */

#endif // RICH_MPI

  class SimData
  {
  public:

    SimData(void):
      width_(read_number("width.txt")),
      init_points_(cartesian_mesh(2*30,2*30,
				  Vector2D(0,0),
				  Vector2D(width_,width_))),
      outer_(0,width_,width_,0),
#ifdef RICH_MPI
      meta_tess_(create_meta_points(outer_),outer_),
      points_(create_mesh_generating_points(60,outer_,meta_tess_)),
      tess_(meta_tess_,points_,outer_),
#else
      tess_(init_points_,outer_),
#endif // RICH_MPI
      pg_(),
      eos_(read_number("adiabatic_index.txt")),
      init_cond_(read_number("ambient_density.txt"),
		 read_number("ambient_pressure.txt"),
		 eos_,
		 read_number("amplitude.txt"),
		 width_),
      rs_(),
      point_motion_(),
      evc_(),
      force_(),
      tsf_(0.3),
      gpg_(),
      sr_(eos_,gpg_),
      fc_(sr_,rs_),
      eu_(),
      cu_(),
      sim_(
#ifdef RICH_MPI
	   meta_tess_,
#endif // RICH_MPI
	   tess_,
	   outer_,
	   pg_,
	   calc_init_cond(tess_,eos_,width_),
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

    double width_;
    vector<Vector2D> init_points_;
    PeriodicBox outer_;
#ifdef RICH_MPI
    VoronoiMesh meta_tess_;
    vector<Vector2D> points_;
#endif // RICH_MPI
    VoronoiMesh tess_;
    SlabSymmetry pg_;
    IdealGas eos_;
    AcousticInitCond init_cond_;
    Hllc rs_;
    Eulerian point_motion_;
    const PeriodicEdgeVelocities evc_;
    ZeroForce force_;
    const SimpleCFL tsf_;
    const PeriodicGhostGenerator gpg_;
    const LinearGaussImproved sr_;
    const ModularFluxCalculator fc_;
    const SimpleExtensiveUpdater eu_;
    const SimpleCellUpdater cu_;
    hdsim sim_;
  };
}

int main(void)
{
#ifdef RICH_MPI
  MPI_Init(NULL, NULL);
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif // RICH_MPI
  try
    {
      SimData sim_data;
      hdsim& sim = sim_data.getSim();
      SafeTimeTermination term_cond(1,1e6);
      WriteTime diag("time.txt");
#ifdef RICH_MPI
      write_snapshot_to_hdf5(sim, "initial_"+int2str(rank)+".h5");
#else
      write_snapshot_to_hdf5(sim, "initial.h5");
#endif // RICH_MPI
      main_loop(sim, 
		term_cond,
		&hdsim::TimeAdvance2Heun,
		&diag);
#ifdef RICH_MPI
      write_snapshot_to_hdf5(sim, "final_"+int2str(rank)+".h5");
#else
      write_snapshot_to_hdf5(sim, "final.h5");
#endif // RICH_MPI
    }
  catch(const UniversalError& eo){
    report_error(eo);
    throw;
  }

#ifdef RICH_MPI
  MPI_Finalize();
#endif // RICH_MPI

  return 0;
}
