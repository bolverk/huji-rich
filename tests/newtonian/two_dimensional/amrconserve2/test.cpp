#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/amr.hpp"
#include "source/newtonian/two_dimensional/modular_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/simple_cfl.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/interpolations/LinearGaussImproved.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/periodic_edge_velocities.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/PeriodicBox.hpp"
#include "source/newtonian/two_dimensional/ghost_point_generators/PeriodicGhostGenerator.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include <boost/random/uniform_int_distribution.hpp>
#include "source/misc/simple_io.hpp"

namespace {
class ConserveRefine : public CellsToRefine
{
public:
	vector<size_t> ToRefine(Tessellation const& tess, vector<ComputationalCell> const& /*cells*/, double time,
	TracerStickerNames const& /*ts*/)const
	{
		boost::random::mt19937_64 gen(static_cast<uint64_t>(time*10000));
		boost::random::uniform_int_distribution<> dist(0, tess.GetPointNo()-1);
		vector<size_t> res(10);
		for (size_t i = 0; i < res.size(); ++i)
		  res[i] = static_cast<size_t>(dist(gen));
		return res;
	}
};

class ConserveRemove : public CellsToRemove
{
public:
	std::pair<vector<size_t>, vector<double> > ToRemove(Tessellation const& tess, vector<ComputationalCell> const& /*cells*/, double time,
	TracerStickerNames const& /*ts*/)const
	{
		boost::random::mt19937_64 gen(static_cast<uint64_t>(time * 11000));
		boost::random::uniform_int_distribution<> dist(0,tess.GetPointNo()-1);
		std::pair<vector<size_t>, vector<double> > res(vector<size_t>(10), vector<double>(10));
		for (size_t i = 0; i < res.first.size(); ++i)
		{
		  res.first[i] = static_cast<size_t>(dist(gen));
			res.second[i] = static_cast<double>(dist(gen));
		}
		return res;
	}
};

double TotalMass(hdsim const& sim)
{
	double res = 0;
	for (int i = 0; i < sim.getTessellation().GetPointNo(); ++i)
		res += sim.getTessellation().GetVolume(i)*sim.getAllCells()[static_cast<size_t>(i)].density;
	return res;
}

vector<ComputationalCell> calc_cells(size_t n)
{
	vector<ComputationalCell> res(n);
	for (size_t i = 0; i < n; ++i)
	{
		res[i].density = 1;
		res[i].pressure = 1;
		res[i].velocity = Vector2D(1, 0);
	}
	return res;
}
}

int main(void)
{
	int np = 30;
	vector<Vector2D> points = cartesian_mesh(np, np, Vector2D(-1, -1), Vector2D(1, 1));
	PeriodicBox outer(Vector2D(-1, -1), Vector2D(1, 1));
	VoronoiMesh tess(points, outer);
	vector<ComputationalCell> cells = calc_cells(static_cast<size_t>(tess.GetPointNo()));
	Hllc rs;
	IdealGas eos(5. / 3.);
	PeriodicGhostGenerator ghost;
	PeriodicEdgeVelocities vedge;
	LinearGaussImproved interp(eos, ghost);
	SimpleCellUpdater cu;
	SimpleExtensiveUpdater eu;
	SimpleCFL tsf(0.3);
	ModularFluxCalculator fc(interp, rs);
	Lagrangian pm;
	SlabSymmetry pg;
	ZeroForce force;

	hdsim sim(tess, outer, pg, cells, eos, pm, vedge, force, tsf, fc, eu, cu);

	ConserveRefine refine;
	ConserveRemove remove;
	ConservativeAMR amr(refine, remove,true, &interp);

	for (size_t i = 0; i < 50; ++i)
	{
		sim.TimeAdvance2Heun();
		amr(sim);
	}
	
	write_number(TotalMass(sim), "mass.txt",9);
	write_snapshot_to_hdf5(sim, "final.h5");

	return 0;
}
