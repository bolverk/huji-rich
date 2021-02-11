#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/interpolations/LinearGaussImproved.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/misc/int2str.hpp"
#include "source/newtonian/two_dimensional/modular_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cfl.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/ColdFlowsExtensiveCalculator.hpp"
#include "source/newtonian/two_dimensional/amr.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"

namespace
{
  vector<ComputationalCell> calc_init_cond(const Tessellation& tess,EquationOfState const& eos)
  {
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    for (size_t i = 0; i < res.size(); ++i)
      {
	res[i].density = 1;
	res[i].pressure = 1e-6;
	Vector2D const& point = tess.GetCellCM(static_cast<int>(i));
	const double r = abs(point);
	res[i].velocity = Vector2D(-point.x / r, -point.y / r);
	//res[i].tracers.push_back(eos.dp2s(1, 1e-6));
	res[i].tracers[0] = eos.dp2s(1, 1e-6);
      }
    return res;
  }

  class NOHGhostGenerator : public GhostPointGenerator
  {
  private:
    EquationOfState const& eos_;
  public:

    NOHGhostGenerator(EquationOfState const& eos) :eos_(eos){}

    boost::container::flat_map<size_t, ComputationalCell> operator() (const Tessellation& tess,
								      const vector<ComputationalCell>& /*cells*/, double time,TracerStickerNames const& /*ts*/) const
    {
      vector<std::pair<size_t, size_t> > outer_edges = GetOuterEdgesIndeces(tess);
      boost::container::flat_map<size_t, ComputationalCell> res;
      for (size_t i = 0; i < outer_edges.size(); ++i)
	{
	  Edge const& edge = tess.GetEdge(static_cast<int>(outer_edges[i].first));
	  size_t ghost_index = static_cast<size_t>(outer_edges[i].second == 1 ? edge.neighbors.first : edge.neighbors.second);
	  ComputationalCell temp;
	  const Vector2D edge_cen = tess.GetMeshPoint(outer_edges[i].second == 0 ? edge.neighbors.first : edge.neighbors.second);
	  const double r = abs(edge_cen);
	  temp.density = (1 + time / r);
	  temp.pressure = 1e-6;
	  temp.velocity = -1.0*edge_cen / r;
	  //temp.tracers.push_back(eos_.dp2s(temp.density, temp.pressure));
	  temp.tracers[0] = eos_.dp2s(temp.density,
				      temp.pressure);
	  res[ghost_index] = temp;
	}
      return res;
    }

    Slope GetGhostGradient(Tessellation const& /*tess*/,
			   vector<ComputationalCell> const& /*cells*/, vector<Slope> const& /*gradients*/,
			   size_t /*ghost_index*/, double /*time*/,Edge const& /*edge*/,TracerStickerNames const& /*ts*/)const
    {
      ComputationalCell temp;
      //temp.tracers.push_back(0);
      temp.tracers[0] = 0;
      return Slope(temp, temp);
    }
  };

  class NohRefine: public CellsToRefine
  {
  private:
    const double maxV_;
		
  public:
    NohRefine(double maxV) :maxV_(maxV){}

    vector<size_t> ToRefine(Tessellation const& tess, vector<ComputationalCell> const& /*cells*/, double /*time*/
			    ,TracerStickerNames const& /*ts*/)const
    {
      vector<size_t> res;
      size_t N = static_cast<size_t>(tess.GetPointNo());
      for (size_t i = 0; i < N; ++i)
	{
	  const double V = tess.GetVolume(static_cast<int>(i));
	  if (V > maxV_)
	    res.push_back(i);
	}
      return res;
    }
  };

  class NohRefineDebug : public CellsToRefine
  {
  public:
    vector<size_t> ToRefine(Tessellation const& /*tess*/, vector<ComputationalCell> const& /*cells*/, double /*time*/
			    ,TracerStickerNames const& /*ts*/)const
    {
      return vector<size_t>();
    }
  };

  class NohRemove : public CellsToRemove
  {
  private:
    const double minV_;
    LinearGaussImproved const& interp_;
  public:
    NohRemove(double minV, LinearGaussImproved const& interp) :minV_(minV), interp_(interp) {}

    std::pair<vector<size_t>, vector<double> > ToRemove(Tessellation const& tess,
							vector<ComputationalCell> const& cells, double /*time*/
							,TracerStickerNames const& /*ts*/)const
    {
      vector<size_t> indeces;
      vector<double> merits;
      size_t N = static_cast<size_t>(tess.GetPointNo());
      for (size_t i = 0; i < N; ++i)
	{
	  const double V = tess.GetVolume(static_cast<int>(i));
	  if (V < minV_)
	    {
	      const double drhox = interp_.GetSlopesUnlimited()[i].xderivative.density;
	      const double drhoy = interp_.GetSlopesUnlimited()[i].yderivative.density;
	      if ((drhox*drhox + drhoy*drhoy)*V < 0.1*cells[i].density*cells[i].density)
		{
		  indeces.push_back(i);
		  merits.push_back(1.0 / V);
		}
	    }
	}
      return std::pair<vector<size_t>, vector<double> >(indeces, merits);
    }
  };


  class NohRemoveDebug : public CellsToRemove
  {
  public:
    std::pair<vector<size_t>, vector<double> > ToRemove(Tessellation const& /*tess*/,
							vector<ComputationalCell> const& /*cells*/, double /*time*/,TracerStickerNames const& /*ts*/)const
    {
      return std::pair<vector<size_t>, vector<double> >();
    }
  };

}

int main(void)
{
  // Set up the initial grid points
  int np = read_int("resolution.txt");
  //int np =30;
  double width = 2;
  // Set up the boundary type for the points
  SquareBox outer(-width*0.5,width*0.5,width*0.5,-width*0.5);

  vector<Vector2D> InitPoints = cartesian_mesh(np, np, Vector2D(-1, -1), Vector2D(1, 1));

  // Set up the tessellation
  VoronoiMesh tess(InitPoints, outer);

  // Set up the Riemann solver
  Hllc rs;

  // Set up the equation of state
  IdealGas eos(read_number("adiabatic_index.txt"));
  //IdealGas eos(5./3.);

  // Set up the point motion scheme
  Lagrangian l_motion;
  RoundCells pointmotion(l_motion, eos,0.5,0.01,true);
  const StationaryBox evc;

  // Set the ghost points
  NOHGhostGenerator ghost(eos);

  // Set up the interpolation
  LinearGaussImproved interpolation(eos, ghost);

  // Set up the external source term
  ZeroForce force;

  ColdFlowsExtensiveCalculator eu(eos,ghost,interpolation);
  SimpleCFL tsf(0.15);
  ModularFluxCalculator fc(interpolation, rs);
  SimpleCellUpdater cu;
  SlabSymmetry pg;

  vector<ComputationalCell> init_cells = calc_init_cond(tess,eos);

  // Set up the simulation
  //#define restart 1

#ifdef restart
  Snapshot snap = read_hdf5_snapshot("c:/sim_data/snap.h5");
  tess.Update(snap.mesh_points);
  init_cells = snap.cells;
#endif
  hdsim sim(tess, outer, pg, init_cells, eos, pointmotion, evc, force, tsf, fc, eu, cu,TracerStickerNames (vector<string> (1,"Entropy"),vector
													   <string>()));
#ifdef restart
  sim.setStartTime(snap.time);
#endif
  // Define the AMR 
  double Vmax = 3 * width*width / (np*np);
  double Vmin = 0.25*width*width / (np*np);
  NohRefine refine(Vmax);
  //NohRefineDebug refine;
  NohRemove remove(Vmin, interpolation);
  //NohRemoveDebug remove;
  ConservativeAMR amr(refine, remove);

  // How long shall we run the simulation?
#ifdef restart
  double tend = 3;
#else
  double tend = 2;
#endif

  // Run main loop of the sim
  while (sim.getTime()<tend)
    {
      try
	{
	  // Advance one time step
	  //			write_snapshot_to_hdf5(sim, "c:/sim_data/snap.h5");
	  sim.TimeAdvance2Heun();
	  write_number(sim.getTime(), "time.txt");
	  amr(sim);
	}
      catch (UniversalError const& eo)
	{
	  DisplayError(eo);
	  throw;
	}
    }

  // Done running the simulation, output the data
  //string filename = "final.h5";
#ifdef restart
  string filename = "c:/sim_data/final2.h5";
#else
  string filename = "final.h5";
#endif
  write_snapshot_to_hdf5(sim, filename);
  return 0;

}

