#ifdef RICH_MPI
#include "source/mpi/MeshPointsMPI.hpp"
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/environment.hpp>
#endif
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "source/misc/int2str.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/tessellation/geometry.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/CustomMotion.hpp"
#include "source/misc/utils.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/two_dimensional/condition_action_sequence.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/amr.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"

using namespace std;

namespace
{

#ifdef RICH_MPI
  vector<Vector2D> process_positions(const SquareBox& boundary)
  {
    const boost::mpi::communicator world;
    const Vector2D lower_left = boundary.getBoundary().first;
    const Vector2D upper_right = boundary.getBoundary().second;
    vector<Vector2D> res(world.size());
    if(world.rank()==0){
      res = RandSquare(world.size(),
		       lower_left.x,upper_right.x,
		       lower_left.y,upper_right.y);
    }
    boost::mpi::broadcast(world, res, 0);
    return res;
  }
#endif // RICH_MPI

	class WedgeNoMove : public CustomMotionCriteria
	{
	public:
		WedgeNoMove(string key, double distanceToWall) :key_(key), distanceToWall_(distanceToWall){}
		bool SatisfyCriteria(size_t index, Tessellation const& tess, vector<ComputationalCell> const& cells,
			double /*time*/)const
		{
			if (cells[index].stickers.at(key_) || tess.GetMeshPoint(static_cast<int>(index)).y
				<(0.1*tess.GetMeshPoint(static_cast<int>(index)).x - 0.428 + distanceToWall_))
				return true;
			else
				return false;
		}
		Vector2D CustomVelocityResult(size_t index, Tessellation const& /*tess*/, vector<ComputationalCell> const& cells,
			double /*time*/)const
		{
			if (!cells[index].stickers.at(key_))
				return cells[index].velocity*0.2;
			else
				return Vector2D(0, 0);
		}
	private:
		string key_;
		double const distanceToWall_;
	};

	vector<ComputationalCell> calc_init_cond
		(const Tessellation& tess)
	{
		vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
		for (size_t i = 0; i < res.size(); ++i){
			res[i].density = 1;
			res[i].pressure = 1;
			res[i].velocity = Vector2D(4, 0);
			const Vector2D r = tess.GetCellCM(static_cast<int>(i));
			res[i].stickers["wedge"] = r.y < 0.1*r.x - 0.428;
		}
		write_number(atan(0.1), "wedge_angle.txt");
		return res;
	}

	class LeftBoundaryEdge : public ConditionActionSequence::Condition
	{
	public:

		LeftBoundaryEdge(void) {}

		pair<bool, bool> operator()
			(const Edge& edge,
			const Tessellation& tess,
			const vector<ComputationalCell>& /*cells*/) const
		{
			if (edge.neighbors.first >= tess.GetPointNo())
			{
				assert(edge.neighbors.second<tess.GetPointNo());
				const Vector2D r = tess.GetMeshPoint(edge.neighbors.second);
				if (r.x>edge.vertices.first.x && r.x > edge.vertices.second.x)
					return pair<bool, bool>(true, false);
				return pair<bool, bool>(false, false);
			}
			if (edge.neighbors.second >= tess.GetPointNo())
			{
				const Vector2D r = tess.GetMeshPoint(edge.neighbors.first);
				if (r.x>edge.vertices.first.x && r.x > edge.vertices.second.x)
					return pair<bool, bool>(true, true);
				return pair<bool, bool>(false, false);
			}
			return pair<bool, bool>(false, false);
		}
	};

	ComputationalCell calc_ghost(void)
	{
		ComputationalCell res;
		res.density = 1;
		res.pressure = 1;
		res.velocity = Vector2D(4, 0);
		res.stickers["wedge"] = false;
		return res;
	}

	Extensive conserved_to_extensive
		(const Conserved& c, const ComputationalCell& cell)
	{
		Extensive res;
		res.mass = c.Mass;
		res.momentum = c.Momentum;
		res.energy = c.Energy;
		for (boost::container::flat_map<string, double>::const_iterator it =
			cell.tracers.begin();
			it != cell.tracers.end(); ++it)
			res.tracers[it->first] = (it->second)*c.Mass;
		return res;
	}

	class ConstantGhost : public ConditionActionSequence::Action
	{
	public:

		ConstantGhost(const ComputationalCell& ghost,
			const RiemannSolver& rs) :
			ghost_(ghost), rs_(rs) {}

		Extensive operator()
		(const Edge& edge,
		 const Tessellation& tess,
		 const Vector2D& /*edge_velocity*/,
		 const vector<ComputationalCell>& cells,
		 const EquationOfState& eos,
		 const bool aux) const
		{
			if (aux)
				assert(edge.neighbors.first < tess.GetPointNo());
			else
				assert(edge.neighbors.second < tess.GetPointNo());
			const Vector2D p = normalize
				(edge.vertices.second -
				edge.vertices.first);
			const Vector2D n = normalize
				(remove_parallel_component
				(aux ?
				edge.vertices.first - tess.GetMeshPoint(edge.neighbors.first) :
				tess.GetMeshPoint(edge.neighbors.second) - edge.vertices.first,
				p));
			const double v = 0;
			const pair<ComputationalCell, ComputationalCell> cc_left_righ =
				aux ?
				pair<ComputationalCell, ComputationalCell>
				(cells.at(static_cast<size_t>(edge.neighbors.first)), ghost_) :
				pair<ComputationalCell, ComputationalCell>
				(ghost_, cells.at(static_cast<size_t>(edge.neighbors.second)));
			const pair<Primitive, Primitive> left_right =
				pair<Primitive, Primitive>
				(convert_to_primitive(cc_left_righ.first, eos),
				convert_to_primitive(cc_left_righ.second, eos));
			const Conserved c = rotate_solve_rotate_back
				(rs_,
				left_right.first,
				left_right.second,
				v, n, p);
			return conserved_to_extensive(c, ghost_);
		}

	private:
		const ComputationalCell ghost_;
		const RiemannSolver& rs_;
	};

	class SimData
	{
	public:

		SimData(double adiabatic_index = 5. / 3.,
			double width = 1) :
			pg_(),
			outer_(-width / 2, width / 2, width / 2, -width / 2),
#ifdef RICH_MPI
			proctess_(process_positions(outer_), outer_),
#endif // RICH_MPI
			tess_(
#ifdef RICH_MPI
			      proctess_,
			      SquareMeshM
			      (30,30,
			       proctess_,
			       outer_.getBoundary().first,
			       outer_.getBoundary().second),
#else
			cartesian_mesh(30, 30,
			outer_.getBoundary().first,
			outer_.getBoundary().second),
#endif
			outer_),
			eos_(adiabatic_index),
			rs_(),
			criteria_("wedge", 0.03),
			lpm_(),
			pm_(lpm_, eos_, outer_),
			point_motion_(pm_, criteria_),
			force_(),
			tsf_(0.3),
			fc_(VectorInitialiser<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*> >
			(pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*>
			(new LeftBoundaryEdge, new ConstantGhost(calc_ghost(), rs_)))
			(pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*>
			(new IsBoundaryEdge, new FreeFlowFlux(rs_)))
			(pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*>
			(new RegularSpecialEdge("wedge"), new RigidWallFlux(rs_)))
			(pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*>
			(new IsBulkEdge, new RegularFlux(rs_)))()),
			eu_(),
			cu_(VectorInitialiser<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> >
			(pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*>
			(new HasSticker("wedge"), new SkipUpdate))()),
			sim_
		  (
#ifdef RICH_MPI
		   proctess_,
#endif // RICH_MPI
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
			cu_)
		{
			write_number(adiabatic_index,
				"adiabatic_index.txt");
			write_number(4.0 / sqrt(5. / 3.),
				"mach_number.txt");
		}

		hdsim& getSim(void)
		{
			return sim_;
		}

	private:
		const SlabSymmetry pg_;
		SquareBox outer_;
#ifdef RICH_MPI
	  VoronoiMesh proctess_;
#endif // RICH_MPI
		VoronoiMesh tess_;
		const IdealGas eos_;
		const Hllc rs_;
		WedgeNoMove criteria_;
		Lagrangian lpm_;
		RoundCells pm_;
		CustomMotion point_motion_;
	  const StationaryBox evc_;
		ZeroForce force_;
		const SimpleCFL tsf_;
		const ConditionActionSequence fc_;
		const SimpleExtensiveUpdater eu_;
		const SimpleCellUpdater cu_;
		hdsim sim_;
	};

	class ObliqueRefine : public CellsToRefine
	{
	private:
		const double maxV_;
		string key_;

	public:
		ObliqueRefine(double maxV,string key) :maxV_(maxV),key_(key) {}

		vector<size_t> ToRefine(Tessellation const& tess, vector<ComputationalCell> const& cells, double /*time*/)const
		{
			vector<size_t> res;
			size_t N = static_cast<size_t>(tess.GetPointNo());
			for (size_t i = 0; i < N; ++i)
			{
				if (tess.GetVolume(static_cast<int>(i)) > maxV_&&!cells[i].stickers.at(key_))
					res.push_back(i);
			}
			return res;
		}
	};
	
	class ObliqueRemove :public CellsToRemove
	{
	private:
		string key_;
		double distanceToWall_;
	public:
		ObliqueRemove(string key, double distanceToWall) :key_(key), distanceToWall_(distanceToWall) {}

		std::pair<vector<size_t>, vector<double> > ToRemove(Tessellation const& tess,
			vector<ComputationalCell> const& cells, double /*time*/)const 
		{
			std::pair<vector<size_t>, vector<double> > res;
			vector<double> merits;
			vector<size_t> indeces;
			size_t N = static_cast<size_t>(tess.GetPointNo());
			for (size_t i = 0; i < N; ++i)
			{
				if (!cells[i].stickers.at(key_))
				{
					Vector2D const& point = tess.GetMeshPoint(static_cast<int>(i));
					if (point.y < (0.1*point.x - 0.428 + distanceToWall_))
					{
						indeces.push_back(i);
						merits.push_back(0.1*point.x - 0.428 + distanceToWall_-point.y);
					}
					else
					{
						if (point.x > (0.5 - distanceToWall_))
						{
							indeces.push_back(i);
							merits.push_back(point.x - (0.5 - distanceToWall_));
						}
					}
				}
			}
			return std::pair<vector<size_t>, vector<double> >(indeces,merits);
		}

	};

	void main_loop(hdsim& sim)
	{
		const int max_iter = 5e6;
		const double tf = 0.5;
		ObliqueRefine refine(0.05*0.05,"wedge");
		ObliqueRemove remove("wedge",0.02);
		NonConservativeAMR amr(refine, remove);

		while (tf > sim.getTime()){
#ifdef RICH_MPI
		  cout << sim.getTime() << endl;
#endif // RICH_MPI
			try{
				sim.TimeAdvance();
				amr(sim);
			}
			catch (UniversalError const& eo){
				DisplayError(eo);
			}

			if (sim.getCycle() > max_iter)
				throw UniversalError
				("Maximum number of iterations exceeded in main loop");
		}
	}
}

int main(void)
{
#ifdef RICH_MPI
  boost::mpi::environment env;
  const boost::mpi::communicator world;
#endif // RICH_MPI

	SimData sim_data;

	hdsim& sim = sim_data.getSim();

	main_loop(sim);

#ifdef RICH_MPI
	write_snapshot_to_hdf5(sim, "process_"+int2str(world.rank())+"_final.h5");
#else
	write_snapshot_to_hdf5(sim, "final.h5");
#endif // RICH_MPI

	return 0;
}


