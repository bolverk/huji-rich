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

  vector<Vector2D> my_convex_hull
  (const Tessellation& tess,
   int index)
  {
    vector<Vector2D> res;
    ConvexHull(res,&tess,index);
    return res;
  }

  vector<Vector2D> distribute_grid
  (const vector<Vector2D>& complete_grid,
   const Tessellation& proc_tess,
   const boost::mpi::communicator& world)
  {
    vector<Vector2D> res;
    const vector<Vector2D> ch_list =
      my_convex_hull(proc_tess,world.get_rank());
    BOOST_FOREACH(const Vector2D& v, complete_grid)
      {
	if(PointInCell(ch_list,v))
	  res.push_back(v);
      }
    return res;
  }
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
      tess_(init_points_,outer_),
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
      sim_(tess_,
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
  try
    {
      SimData sim_data;
      hdsim& sim = sim_data.getSim();
      SafeTimeTermination term_cond(1,1e6);
      WriteTime diag("time.txt");
      write_snapshot_to_hdf5(sim, "initial.h5");
      main_loop(sim, 
		term_cond,
		&hdsim::TimeAdvance2Heun,
		&diag);
      write_snapshot_to_hdf5(sim, "final.h5");
    }
  catch(const UniversalError& eo){
    report_error(eo);
    throw;
  }

  return 0;
}
