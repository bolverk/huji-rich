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
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
//#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/source_terms/ConservativeForce.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/misc/mesh_generator.hpp"

using namespace std;
using namespace simulation2d;

namespace {
class PeriodicDriver: public Acceleration
{
public:

  PeriodicDriver(double wavelenth,
		 double amplitude,
		 double phase_velocity):
    k_(2*M_PI/wavelenth),
    amp_(amplitude),
    v_(phase_velocity) {}

  Vector2D Calculate
  (Tessellation const& tess,
   vector<Primitive> const& /*cells*/,
   int point,
   vector<Conserved> const& /*fluxes*/,
   vector<Vector2D> const& /*point_velocity*/,
   HydroBoundaryConditions const& /*hbc*/,
   vector<vector<double> > const& /*tracers*/,
   double t,
   double /*dt*/)
  {
    const double x = tess.GetMeshPoint(point).x;
    const double acceleration = amp_*sin(k_*x)*sin(k_*v_*t);
    return Vector2D(-acceleration,0);
  }

private:
  const double k_;
  const double amp_;
  const double v_;
};

  namespace {
    vector<ComputationalCell> calc_init_cond(const Tessellation& tess)
    {
      const double density = read_number("mean_density.txt");
      const double pressure = read_number("mean_pressure.txt");
      vector<ComputationalCell> res(tess.GetPointNo());
      for(size_t i=0;i<res.size();++i){
	res[i].density = density;
	res[i].pressure = pressure;
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
    init_points_(cartesian_mesh(30,30,Vector2D(0,0),
				Vector2D(width_,width_))),			       
    outer_(0,width_,width_,0),
    pg_(),
    tess_(init_points_, outer_),
    eos_(read_number("adiabatic_index.txt")),
    pm_naive_(),
    rs_(),
    hbc_(rs_),
    //    point_motion_(pm_naive_,hbc_),
    acc_(read_number("wavelength.txt"),
	 read_number("amplitude.txt"),
	 read_number("phase_velocity.txt")),
	force_(acc_),
    sim_(tess_,
	 outer_,
	 pg_,
	 calc_init_cond(tess_),
	 eos_,
	 pm_naive_,
	 force_,
	 hbc_) {}

  hdsim& getSim(void)
  {
    return sim_;
  }

private:
  const double width_;
  const vector<Vector2D> init_points_;
  const SquareBox outer_;
  const SlabSymmetry pg_;
  VoronoiMesh tess_;
  const IdealGas eos_;
  Lagrangian pm_naive_;
  const Hllc rs_;
  const RigidWallHydro hbc_;
  //  RoundCells point_motion_;
  PeriodicDriver acc_;
  ConservativeForce force_;
  hdsim sim_;
};

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(1, 1e6);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      &hdsim::TimeAdvance,
	      &diag);
  }
}

int main(void)
{
  SimData sim_data;
  hdsim& sim = sim_data.getSim();

  my_main_loop(sim);
  
  write_snapshot_to_hdf5(sim, "final.h5");

  return 0;
}

