#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include "source/tessellation/geometry.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/interpolations/pcm2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/step2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/test_2d/random_pert.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/misc/mesh_generator.hpp"

using namespace std;
using namespace simulation2d;

namespace {
  class SimData
  {
  public:

    SimData(void):
      width_(1),
      np_(10),
      init_points_(RandSquare(np_*np_,0,width_,0,width_)),
      outer_(0, width_, width_, 0),
      tess_(),
      interpm_(),
      density_(1),
      pressure_(1),
      xvelocity_(0),
      yvelocity_(0),
      eos_(5./3.),
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

  private:

    const double width_;
    const int np_;
    const vector<Vector2D> init_points_;
    const SquareBox outer_;
    VoronoiMesh tess_;
    PCM2D interpm_;
    const Uniform2D density_;
    const Uniform2D pressure_;
    const Uniform2D xvelocity_;
    const Uniform2D yvelocity_;
    const IdealGas eos_;
    Eulerian bpm_;
    const Hllc rs_;
    const RigidWallHydro hbc_;
    RoundCells point_motion_;
    ZeroForce force_;
    hdsim sim_;
  };

  vector<double> volume_list(hdsim const& sim)
  {
    vector<double> res(sim.GetCellNo(),0);
    for(int i=0;i<sim.GetCellNo();++i)
      res[i] = sim.GetCellVolume(i);
    return res;
  }

  class WriteMinMaxVolume: public DiagnosticFunction
  {
  public:

    WriteMinMaxVolume(string const& fname):
      fhandle_(fname.c_str()) {}

    void operator()(hdsim const& sim)
    {
      fhandle_ << min(volume_list(sim)) << " "
	       << max(volume_list(sim)) << endl;
    }

    ~WriteMinMaxVolume(void)
    {
      fhandle_.close();
    }

  private:

    ofstream fhandle_;
  };

  void my_main_loop(hdsim& sim)
  {
    CycleTermination term_cond(100);
    WriteMinMaxVolume diag("res.txt");
    main_loop(sim,
	      term_cond,
	      &hdsim::TimeAdvance,
	      &diag);
  }
}

int main(void)
{
  try{
    SimData sim_data;
    hdsim& sim = sim_data.getSim();

    my_main_loop(sim);
  }
  catch(UniversalError const& eo){
    DisplayError(eo);
  }

  return 0;
}

