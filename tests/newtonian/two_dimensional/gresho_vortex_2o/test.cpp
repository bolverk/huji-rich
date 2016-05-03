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
#include "source/newtonian/two_dimensional/modular_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/ghost_point_generators/RigidWallGenerator.hpp"
#include "source/newtonian/two_dimensional/interpolations/LinearGaussImproved.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"

using namespace std;
using namespace simulation2d;

namespace {

  double calc_vq2r(double r)
  {
    if(r<0.2)
      return 5;
    else if(r>0.4)
      return 0;
    return 2.0/r-5.0;
  }

  double azimuthal_velocity(double r)
  {
    if(r<0.2)
      return 5*r;
    else if(r>0.4)
      return 0;
    else
      return 2-5*r;
  }

  double calc_pressure(double r)
  {
    if(r<0.2)
      return 5+(25./2.)*pow(r,2);
    else if(r>0.4)
      return 3+4*log(2.);
    else
      return 9+(25./2.)*pow(r,2)-20*r+4*log(r/0.2);
  }

  vector<ComputationalCell> calc_init_cond(const Tessellation& tess)
  {
    vector<ComputationalCell> res
      (static_cast<size_t>(tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i){
      const Vector2D r = tess.GetCellCM(static_cast<int>(i));
      const double radius = abs(r);
      res[i].density = 1;
      res[i].pressure = calc_pressure(radius);
      res[i].velocity = calc_vq2r(radius)*Vector2D(-r.y,r.x);
    }
    return res;
  }

  class VelocityX: public SpatialDistribution
  {
  public:
    double operator()(Vector2D const& r) const
    {
      if(abs(r)<-0)
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
int GetWS(void)
{
	int ws;
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	return ws;
}

  vector<Vector2D> process_positions(const SquareBox& boundary)
  {
	  int ws,rank;
	  MPI_Comm_size(MPI_COMM_WORLD, &ws);
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const Vector2D lower_left = boundary.getBoundary().first;
    const Vector2D upper_right = boundary.getBoundary().second;
    vector<Vector2D> res(ws);
	vector<double> temp(static_cast<size_t>(ws) * 2);
    if(rank==0)
	{
      res = RandSquare(ws,
		       lower_left.x,upper_right.x,
		       lower_left.y,upper_right.y);
	  temp = list_serialize(res);
    }
	MPI_Bcast(&temp[0], ws * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	Vector2D vtemp;
	res = list_unserialize(temp, vtemp);
    return res;
  }
#endif

class SimData
{
public:
  
  SimData(void):
#ifdef RICH_MPI
	  ws_(GetWS()),
#endif
    outer_(-0.5,0.5,0.5,-0.5),
    rs_(),
#ifdef RICH_MPI
	  proctess_(process_positions(outer_), outer_),
#endif
    tess_(
#ifdef RICH_MPI
		proctess_,
		SquareMeshM(30*ws_,30*ws_,proctess_,outer_.getBoundary().first,
			outer_.getBoundary().second),
#else
		cartesian_mesh
		(30, 30,
			outer_.getBoundary().first,
			outer_.getBoundary().second),
#endif
	  outer_),
    eos_(5./3.),
    bpm_(),
    point_motion_(bpm_,eos_),
    evc_(),
    force_(),
    tsf_(0.3),
    gpg_(),
    sr_(eos_,gpg_),
    fc_(sr_,rs_),
    eu_(),
    sim_(
#ifdef RICH_MPI
		proctess_,
#endif
		tess_,
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
#ifdef RICH_MPI
	const int ws_;
#endif
  const SlabSymmetry pg_;
  const SquareBox outer_;
  const Hllc rs_;
#ifdef RICH_MPI
  VoronoiMesh proctess_;
#endif
  VoronoiMesh tess_;
  const IdealGas eos_;
  Lagrangian bpm_;
  RoundCells point_motion_;
  const StationaryBox evc_;
  ZeroForce force_;
  const SimpleCFL tsf_;
  const RigidWallGenerator gpg_;
  const LinearGaussImproved sr_;
  const ModularFluxCalculator fc_;
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
	      &hdsim::TimeAdvance2Heun,
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
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  write_snapshot_to_hdf5(sim, "process_"+int2str(rank)+"_final.h5");
  MPI_Finalize();
  #else
  write_snapshot_to_hdf5(sim, "final.h5");
  #endif

  return 0;
}
