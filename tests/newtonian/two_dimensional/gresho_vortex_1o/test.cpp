#ifdef RICH_MPI
#include "source/mpi/MeshPointsMPI.hpp"
#endif
#include <iostream>
#include <cmath>
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"

using namespace std;
using namespace simulation2d;

namespace {

  /*
  double calc_vq2r(double r)
  {
    if(r<0.2)
      return 5;
    else if(r>0.4)
      return 0;
    return 2.0/r-5.0;
  }
  */
  
  double azimuthal_velocity(double r)
  {
    if(r<0.2)
      return 5*r;
    else if(r>0.4)
      return 0;
    else
      return 2-5*r;
  }

  class VelocityX: public SpatialDistribution
  {
  public:
    double operator()(Vector2D const& r) const
    {
      if(abs(r)<=0)
	return 0;
      else
	return -azimuthal_velocity(abs(r))*r.y/abs(r);
    }
  };

  class VelocityY: public SpatialDistribution
  {
  public:
    double operator()(Vector2D const& r) const
    {
      if(abs(r)<=0)
	return 0;
      else
	return azimuthal_velocity(abs(r))*r.x/abs(r);
    }
  };

  double calc_pressure(double r)
  {
    if(r<0.2)
      return 5+(25./2.)*pow(r,2);
    else if(r>0.4)
      return 3+4*log(2.);
    else
      return 9+(25./2.)*pow(r,2)-20*r+4*log(r/0.2);
  }

  /*
  vector<ComputationalCell> calc_init_cond
  (const Tessellation& tess)
  {
    vector<ComputationalCell> res
      (static_cast<size_t>
       (tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i){
      const Vector2D r = tess.GetCellCM
	(static_cast<int>(i));
      const double radius = abs(r);
      res[i].density = 1;
      res[i].pressure = calc_pressure(radius);
      res[i].velocity = calc_vq2r(radius)*Vector2D(-r.y,r.x);
    }
    return res;
  }
  */

  class Pressure: public SpatialDistribution
  {
  public:
    double operator()(Vector2D const& r) const
    {
      return calc_pressure(abs(r));
    }
  };
}

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
  vector<ComputationalCell> calc_init_cond
  (const Tessellation& tess)
  {
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i){
      res[i].density = 1;
      const Vector2D r_vec = tess.GetMeshPoint(static_cast<int>(i));
      const double r = abs(r_vec);
      res[i].pressure = calc_pressure(r);
      res[i].velocity = -azimuthal_velocity(r)*zcross(r_vec)/abs(r);
    }
    return res;
  } 
}

class SimData
{
public:
  
  SimData(void):
    pg_(),
    outer_(-0.5,0.5,0.5,-0.5),
    rs_(),
    tess_(cartesian_mesh(30,30,outer_.getBoundary().first,
			 outer_.getBoundary().second),
	  outer_),
    eos_(5./3.),
    bpm_(),
    point_motion_(bpm_,eos_),
    evc_(),
    force_(),
    tsf_(0.3),
    fc_(rs_),
    eu_(),
    cu_(),
    sim_(tess_,
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
  const SlabSymmetry pg_;
  const SquareBox outer_;
  #ifdef RICH_MPI
  VoronoiMesh proc_tess_;
  #endif
  const Hllc rs_;
  VoronoiMesh tess_;
  const IdealGas eos_;
  Lagrangian bpm_;
  RoundCells point_motion_;
  const StationaryBox evc_;
  const ZeroForce force_;
  const SimpleCFL tsf_;
  const SimpleFluxCalculator fc_;
  const SimpleExtensiveUpdater eu_;
  const SimpleCellUpdater cu_;
  hdsim sim_;
};

namespace {

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.003, 1e6);
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

  write_snapshot_to_hdf5(sim, "initial.h5");

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
