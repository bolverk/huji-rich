#include <iostream>
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/tessellation/shape_2d.hpp"
#include "source/newtonian/test_2d/piecewise.hpp"
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"

using namespace std;
using namespace simulation2d;

namespace {

  vector<ComputationalCell> calc_init_cond(const Tessellation& tess)
  {
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i){
      res[i].density = 1;
      res[i].pressure = abs(tess.GetMeshPoint(static_cast<int>(i)))<0.05 ?
	1e5 : 1;
      res[i].velocity = Vector2D(0,0);
    }
    return res;
  }

  class SimData
  {
  public:
    SimData(Vector2D lower_left = Vector2D(0.02,0)+Vector2D(0,0),
	    Vector2D upper_right = Vector2D(0.02,0)+Vector2D(1,1)):
      pg_(Vector2D(0,0), Vector2D(0,1)),
      mesh_(cartesian_mesh(50,50,lower_left,upper_right)),
      outer_(lower_left, upper_right),
      tess_(mesh_,outer_),
      eos_(5./3.),
      rs_(),
      raw_point_motion_(),
      point_motion_(raw_point_motion_,eos_),
      evc_(),
      force_(pg_.getAxis()),
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
    const CylindricalSymmetry pg_;
    const vector<Vector2D> mesh_;
    const SquareBox outer_;
    VoronoiMesh tess_;
    const IdealGas eos_;
    const Hllc rs_;
    Lagrangian raw_point_motion_;
    RoundCells point_motion_;
    const StationaryBox evc_;
    CylindricalComplementary force_;
    const SimpleCFL tsf_;
    const SimpleFluxCalculator fc_;
    const SimpleExtensiveUpdater eu_;
    const SimpleCellUpdater cu_;
    hdsim sim_;
  };

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.01,1e6);
    WriteTime diag("time.txt");
    //TotalConservedHistory diag("history.txt");
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
