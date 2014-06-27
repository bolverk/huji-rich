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
#include "source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/step2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
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

using namespace std;
using namespace simulation2d;

/*! \brief Contains all the data the simulation needs
 */
class SimData 
{
public:

  SimData(void):
    width_(1),
    //    init_points_(square_grid(width_,30)),
    init_points_(cartesian_mesh(30,30,Vector2D(0,0),
				Vector2D(width_,width_))),
    outer_(0,width_,width_,0),
    tess_(),
    eos_(read_number("adiabatic_index.txt")),
    interpm_(eos_,outer_,hbc_,true,false),
    density_1d_(1,10,0.3,0.5),
    density_(density_1d_),
    pressure_1d_(density_1d_,1,
		 eos_.getAdiabaticIndex()),
    pressure_(pressure_1d_),
    sound_speed_(pressure_1d_,density_1d_,eos_),
    c_edge_(sound_speed_.EvalAt(width_)),
    rv_edge_(calc_riemann_invariant
	     (0,c_edge_,
	      eos_.getAdiabaticIndex(),0)),
    xvelocity_1d_(rv_edge_,0,sound_speed_,
		  eos_.getAdiabaticIndex()),
    xvelocity_(xvelocity_1d_),
    yvelocity_(0),
    bpm_(),
    rs_(),
    hbc_(rs_),
    point_motion_(bpm_,hbc_),
    force_(),
    sim_(init_points_,
	 tess_,
	 interpm_,
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

  ~SimData(void) {}

private:

  const double width_;
  const vector<Vector2D> init_points_;
  const SquareBox outer_;
  VoronoiMesh tess_;
  const IdealGas eos_;
  LinearGaussConsistent interpm_;
  const Collela density_1d_;
  const Profile1D density_;
  const ConstEntropy pressure_1d_;
  const Profile1D pressure_;
  const SoundSpeedDist sound_speed_;
  const double c_edge_;
  const double rv_edge_;
  const ConstRiemannInv xvelocity_1d_;
  const Profile1D xvelocity_;
  const Uniform2D yvelocity_;
  Lagrangian bpm_;
  const Hllc rs_;
  const RigidWallHydro hbc_;
  RoundCells point_motion_;
  ZeroForce force_;
  hdsim sim_;
};

namespace {

  void my_main_loop(hdsim& sim)
  {
    SafeTimeTermination term_cond(0.02, 1e6);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      2,
	      &diag);
  }
}

namespace {
void report_error(UniversalError const& eo)
{
cout << eo.GetErrorMessage() << endl;
cout.precision(14);
for(int i=0;i<(int)eo.GetFields().size();++i)
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

