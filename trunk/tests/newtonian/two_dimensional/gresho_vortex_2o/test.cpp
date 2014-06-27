#ifdef RICH_MPI
#include "source/mpi/MeshPointsMPI.hpp"
#endif
#include <iostream>
#include <cmath>
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
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

using namespace std;
using namespace simulation2d;

namespace {
  
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
      if(abs(r)==0)
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
      if(abs(r)==0)
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

class SimData
{
public:
  
  SimData(void):
    outer_(-0.5,0.5,0.5,-0.5),
    #ifdef RICH_MPI
    proc_tess_(process_positions(outer_),outer_),
    #endif
    rs_(),
    hbc_(rs_),
    tess_(),
    eos_(5./3.),
    interpm_(eos_,outer_,hbc_,true,false),
    bpm_(),
    point_motion_(bpm_,hbc_),
    force_(),
    sim_(
	 #ifdef RICH_MPI
	 distribute_grid(proc_tess_,
			 CartesianGridGenerator
			 (30,30,outer_.getBoundary().first,
			  outer_.getBoundary().second)),
	 #else
	 cartesian_mesh(30,30,outer_.getBoundary().first,
			outer_.getBoundary().second),
	 #endif
	 tess_,
	 #ifdef RICH_MPI
	 proc_tess_,
	 #endif
	 interpm_,
	 Uniform2D(1),
	 Pressure(),
	 VelocityX(),
	 VelocityY(),
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
  
  const SquareBox outer_;
  #ifdef RICH_MPI
  VoronoiMesh proc_tess_;
  #endif
  const Hllc rs_;
  const RigidWallHydro hbc_;
  VoronoiMesh tess_;
  const IdealGas eos_;
  LinearGaussConsistent interpm_;
  Lagrangian bpm_;
  RoundCells point_motion_;
  ZeroForce force_;
  hdsim sim_;
};

namespace {

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.003, 1e6);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      2,
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
