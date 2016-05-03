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
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/test_1d/simple_waves_ideal_gas.hpp"
#include "source/newtonian/test_2d/profile_1d.hpp"
#include "source/newtonian/test_1d/collela.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/two_dimensional/ghost_point_generators/RigidWallGenerator.hpp"
#include "source/newtonian/two_dimensional/interpolations/LinearGaussImproved.hpp"
#include "source/newtonian/two_dimensional/modular_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"

using namespace std;
using namespace simulation2d;

namespace {
  vector<ComputationalCell> calc_init_cond
  (const Tessellation& tess,
   const IdealGas& eos)
  {
    vector<ComputationalCell> res
      (static_cast<size_t>(tess.GetPointNo()));
    const Collela density_prof(1,10,0.3,0.5);
    const ConstEntropy pressure_prof
      (density_prof,
       1, eos.getAdiabaticIndex());
    const SoundSpeedDist sound_speed_prof
      (pressure_prof, density_prof, eos);
    const double c_edge =
      sound_speed_prof(1);
    const double rv_edge =
      calc_riemann_invariant
      (0,c_edge,eos.getAdiabaticIndex(),0);
    const ConstRiemannInv xvelocity_prof
      (rv_edge,0,sound_speed_prof,
       eos.getAdiabaticIndex());
    for(size_t i=0;i<res.size();++i){
      const Vector2D r =
	tess.GetCellCM(static_cast<int>(i));
      res[i].density = density_prof(r.x);
      res[i].pressure = pressure_prof(r.x);
      res[i].velocity = Vector2D
	(xvelocity_prof(r.x),0);
    }
    return res;
  }
}

/*! \brief Contains all the data the simulation needs
 */
class SimData 
{
public:

  SimData(void):
    width_(1),
    outer_(0,width_,width_,0),
    tess_(cartesian_mesh(30,30,Vector2D(0,0),
			 Vector2D(width_,width_)),
	  outer_),
    eos_(read_number("adiabatic_index.txt")),
    bpm_(),
    rs_(),
    gpg_(),
    sr_(eos_,gpg_),
    point_motion_(bpm_,eos_),
    evc_(),
    force_(),
    tsf_(0.3),
    fc_(sr_,rs_),
    eu_(),
    sim_(tess_,
	 outer_,
	 pg_,
	 calc_init_cond(tess_,eos_),
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

  ~SimData(void) {}

private:

  const double width_;
  const SquareBox outer_;
  const SlabSymmetry pg_;
  VoronoiMesh tess_;
  const IdealGas eos_;
  Lagrangian bpm_;
  const Hllc rs_;
  const RigidWallGenerator gpg_;
  const LinearGaussImproved sr_;
  RoundCells point_motion_;
  const StationaryBox evc_;
  ZeroForce force_;
  const SimpleCFL tsf_;
  const ModularFluxCalculator fc_;
  const SimpleExtensiveUpdater eu_;
  const SimpleCellUpdater cu_;
  hdsim sim_;
};

namespace {

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.02, 1e6);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      &hdsim::TimeAdvance2Heun,
	      &diag);
  }
}

namespace {
void report_error(UniversalError const& eo)
{
cout << eo.GetErrorMessage() << endl;
cout.precision(14);
for(size_t i=0;i<eo.GetFields().size();++i)
  cout << eo.GetFields()[i] << " = "
       << eo.GetValues()[i] << endl;
}
}

int main(void)
{
  SimData sim_data;
  hdsim& sim = sim_data.getSim();


  write_snapshot_to_hdf5(sim, "initial.h5");

try{
my_main_loop(sim);
}
 catch(UniversalError const& eo)
   {
report_error(eo);
throw;
}

  write_snapshot_to_hdf5(sim, "final.h5");

  return 0;
}

