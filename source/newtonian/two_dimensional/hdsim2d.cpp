#include <cmath>
#include <algorithm>
#include "hdsim2d.hpp"
#include "hydrodynamics_2d.hpp"
#ifdef RICH_MPI
#include <boost/mpi/communicator.hpp>
#endif //RICH_MPI

using namespace std;

namespace
{
	class CellEdgesGetter : public LazyList<Edge>
	{
	public:

		CellEdgesGetter(const Tessellation& tess, int n) :
			tess_(tess), edge_indices_(tess.GetCellEdges(n)) {}

		size_t size(void) const
		{
			return edge_indices_.size();
		}

		Edge operator[](size_t i) const
		{
			return tess_.GetEdge(edge_indices_[i]);
		}

	private:
		const Tessellation& tess_;
		const vector<int> edge_indices_;
	};
}

namespace
{

	vector<Extensive> init_extensives(const Tessellation& tess,
		const PhysicalGeometry& pg,
		const vector<ComputationalCell>& cells,
		const EquationOfState& eos)
	{
		vector<Extensive> res(cells.size());
		for (size_t i = 0; i<cells.size(); ++i){
			const ComputationalCell& cell = cells[i];
			const double volume =
				pg.calcVolume
				(serial_generate(CellEdgesGetter(tess, static_cast<int>(i))));
			const double mass = volume*cell.density;
			res[i].mass = mass;
			res[i].energy = eos.dp2e(cell.density, cell.pressure, cell.tracers)*mass +
				0.5*mass*ScalarProd(cell.velocity, cell.velocity);
			res[i].momentum = mass*cells[i].velocity;
			for (boost::container::flat_map<std::string, double>::const_iterator it =
				cells[i].tracers.begin();
				it != cells[i].tracers.end(); ++it)
				res[i].tracers[it->first] = (it->second)*mass;
		}
		return res;
	}
}

hdsim::hdsim
(
#ifdef RICH_MPI
	Tessellation& proctess,
#endif
	Tessellation& tess,
 const OuterBoundary& obc,
 const PhysicalGeometry& pg,
 const vector<ComputationalCell>& cells,
 const EquationOfState& eos,
 const PointMotion& point_motion,
 const EdgeVelocityCalculator& evc,
 const SourceTerm& source,
 const TimeStepFunction& tsf,
 const FluxCalculator& fc,
 const ExtensiveUpdater& eu,
 const CellUpdater& cu) :
#ifdef RICH_MPI
	proctess_(proctess),
#endif
  tess_(tess),
  obc_(obc),
  eos_(eos),
  cells_(cells),
  extensives_
  (init_extensives
   (tess,
    pg,
    cells,
    eos)),
  point_motion_(point_motion),
  edge_velocity_calculator_(evc),
  source_(source),
  time_(0),
  cycle_(0),
  pg_(pg),
  tsf_(tsf),
  fc_(fc),
  eu_(eu),
  cu_(cu),
  cache_data_(tess, pg) {}

hdsim::~hdsim(void) {}

namespace {

#ifdef RICH_MPI
	void exchange_cells(const Tessellation& tess, vector<ComputationalCell>& cells)
	{
		const boost::mpi::communicator world;
		const vector<int>& correspondents = tess.GetDuplicatedProcs();
		const vector<vector<int> >& duplicated_points = tess.GetDuplicatedPoints();
		vector<vector<ComputationalCell> > incoming(correspondents.size());
		vector<boost::mpi::request> requests;
		for (size_t i = 0; i < correspondents.size(); ++i)
		{
			requests.push_back(world.isend(correspondents[i], 0, VectorValues(cells, duplicated_points[i])));
			requests.push_back(world.irecv(correspondents[i], 0, incoming[i]));
		}
		const vector<vector<int> >& ghost_indices = tess.GetGhostIndeces();
		boost::mpi::wait_all(requests.begin(), requests.end());
		cells.resize(tess.GetTotalPointNumber());
		for (size_t i = 0; i < incoming.size(); ++i)
		{
			for (size_t j = 0; j < incoming.at(i).size(); ++j)
			{
				cells.at(ghost_indices.at(i).at(j)) = incoming.at(i).at(j);
			}
		}
	}

	void exchange_point_velocities
		(const Tessellation& tess,
			vector<Vector2D>& point_velocities)
	{
		const boost::mpi::communicator world;
		const vector<int>& correspondents = tess.GetDuplicatedProcs();
		const vector<vector<int> >& duplicated_points =
			tess.GetDuplicatedPoints();
		vector<vector<Vector2D> > incoming(correspondents.size());
		vector<boost::mpi::request> requests;
		for (size_t i = 0; i < correspondents.size(); ++i) {
			requests.push_back
				(world.isend
					(correspondents[i], 0, VectorValues
						(point_velocities, duplicated_points[i])));
			requests.push_back
				(world.irecv
					(correspondents[i], 0, incoming[i]));
		}
		const vector<vector<int> >& ghost_indices =
			tess.GetGhostIndeces();
		boost::mpi::wait_all(requests.begin(), requests.end());
		point_velocities.resize(tess.GetTotalPointNumber());
		for (size_t i = 0; i < incoming.size(); ++i) {
			for (size_t j = 0; j < incoming.at(i).size(); ++j) {
				point_velocities.at(ghost_indices.at(i).at(j)) =
					incoming.at(i).at(j);
			}
		}
	}
#endif // RICH_MPI
}

void hdsim::TimeAdvance(void)
{
	vector<Vector2D> point_velocities =
		point_motion_(tess_, cells_, time_);

#ifdef RICH_MPI
	exchange_point_velocities(tess_, point_velocities);
	exchange_cells(tess_, cells_);
#endif

	vector<Vector2D> edge_velocities =
	  edge_velocity_calculator_(tess_,point_velocities);

	const double dt = tsf_(tess_,
		cells_,
		eos_,
		edge_velocities,
		time_);

	point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities);
	edge_velocities =
	  edge_velocity_calculator_(tess_,point_velocities);

	const vector<Extensive> fluxes =
	  fc_
	  (tess_,
	   edge_velocities,
	   cells_,
	   extensives_,
	   cache_data_,
	   eos_,
	   time_,
	   dt);


	//  update_extensives(fluxes,
	eu_(fluxes,
		pg_,
		tess_,
		dt,
		cache_data_,
		cells_,
		extensives_);

	ExternalForceContribution(tess_,
		pg_,
		cache_data_,
		cells_,
		fluxes,
		point_velocities,
		source_,
		time_,
		dt,
		extensives_);

#ifdef RICH_MPI
	MoveMeshPoints(point_velocities, dt, tess_,proctess_);
#else
	MoveMeshPoints(point_velocities, dt, tess_);
#endif
	cache_data_.reset();

	cells_ = cu_(tess_, pg_, eos_, extensives_, cells_,
		cache_data_);

	time_ += dt;
	cycle_++;
}

namespace 
{
	vector<Extensive> average_extensive
		(const vector<Extensive>& extensives_1,
		const vector<Extensive>& extensives_2)
	{
		assert(extensives_1.size() == extensives_2.size());
		vector<Extensive> res(extensives_1.size());
		for (size_t i = 0; i<extensives_1.size(); ++i)
			res[i] = 0.5*(extensives_1[i] + extensives_2[i]);
		return res;
	}
}

void hdsim::TimeAdvance2Heun(void)
{
	vector<Vector2D> point_velocities = point_motion_(tess_, cells_, time_);

	vector<Vector2D> edge_velocities =
	  edge_velocity_calculator_(tess_,point_velocities);

	const double dt = tsf_(tess_, cells_, eos_, edge_velocities, time_);

	point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities);
	edge_velocities =
	  edge_velocity_calculator_(tess_,point_velocities);

	const vector<Extensive> mid_fluxes =
	  fc_
	  (tess_,
	   edge_velocities,
	   cells_,
	   extensives_,
	   cache_data_,
	   eos_,
	   time_,
	   dt);

	vector<Extensive> mid_extensives = extensives_;
	eu_(mid_fluxes, pg_, tess_, dt, cache_data_, cells_, mid_extensives);

	ExternalForceContribution
	  (tess_,
	   pg_,
	   cache_data_,
	   cells_,
	   mid_fluxes,
	   point_velocities,
	   source_,
	   time_,
	   dt,
	   mid_extensives);

	MoveMeshPoints(point_velocities, dt, tess_);
	cache_data_.reset();

	const vector<ComputationalCell> mid_cells = cu_(tess_, pg_, eos_, mid_extensives, cells_, cache_data_);

	edge_velocities =
	  edge_velocity_calculator_(tess_,point_velocities);

	const vector<Extensive> fluxes =
	  fc_
	  (tess_,
	   edge_velocities,
	   mid_cells,
	   mid_extensives,
	   cache_data_,
	   eos_,
	   time_,
	   dt);

	eu_(fluxes, pg_, tess_, dt, cache_data_, cells_, extensives_);

	ExternalForceContribution
	  (tess_,
	   pg_,
	   cache_data_,
	   mid_cells,
	   fluxes,
	   point_velocities,
	   source_,
	   time_+dt,
	   dt,
	   extensives_);

	extensives_ = average_extensive(extensives_, mid_extensives);

	cells_ = cu_(tess_, pg_, eos_, extensives_, cells_, cache_data_);

	time_ += dt;
	++cycle_;
}

namespace {

	template<class T> class AverageCalculator : public LazyList<T>
	{
	public:

		AverageCalculator(const vector<T>& ll1,
			const vector<T>& ll2) :
			ll1_(ll1), ll2_(ll2)
		{
			assert(ll1.size() == ll2.size());
		}

		size_t getLength(void) const
		{
			return ll1_.size();
		}

		T operator()(size_t i) const
		{
			return 0.5*(ll1_[i] + ll2_[i]);
		}

	private:
		const vector<T>& ll1_;
		const vector<T>& ll2_;
	};

	template<class T> vector<T> average(const vector<T>& v1,
		const vector<T>& v2)
	{
		return serial_generate(AverageCalculator<T>(v1, v2));
	}
}

const PhysicalGeometry& hdsim::getPhysicalGeometry(void) const
{
	return pg_;
}

Tessellation& hdsim::getTessellation(void) 
{
	return tess_;
}

const Tessellation& hdsim::getTessellation(void) const
{
	return tess_;
}


void hdsim::addTracer(const string& name,
	const SpatialDistribution& tp)
{
	for (size_t i = 0; i<cells_.size(); ++i)
		cells_[i].tracers[name] = tp(tess_.GetMeshPoint(static_cast<int>(i)));
}

// Diagnostics

double hdsim::getTime(void) const
{
	return time_;
}

int hdsim::getCycle(void) const
{
	return cycle_;
}

const vector<ComputationalCell>& hdsim::getAllCells(void) const
{
	return cells_;
}

vector<ComputationalCell>& hdsim::getAllCells(void)
{
	return cells_;
}

void hdsim::recalculatePrimitives(void)
{
	cells_ = cu_(tess_, pg_, eos_, extensives_, cells_,
		cache_data_);
}

void hdsim::recalculateExtensives(void)
{
	for (size_t i = 0; i<extensives_.size(); ++i){
		const ComputationalCell& cell = cells_[i];
		const double volume = cache_data_.volumes[i];
		const double mass = volume*cell.density;
		extensives_[i].mass = mass;
		extensives_[i].energy = eos_.dp2e(cell.density, cell.pressure, cell.tracers)*mass +
			0.5*mass*ScalarProd(cell.velocity, cell.velocity);
		extensives_[i].momentum = mass*cell.velocity;
		for (boost::container::flat_map<std::string, double>::const_iterator it =
			cell.tracers.begin();
			it != cell.tracers.end(); ++it)
			extensives_[i].tracers[it->first] = (it->second)*mass;
	}
}
void hdsim::setStartTime(double t_start)
{
	time_ = t_start;
}

const EquationOfState& hdsim::getEos(void) const
{
	return eos_;
}

const OuterBoundary& hdsim::getOuterBoundary(void) const
{
	return obc_;
}

const vector<Extensive>& hdsim::getAllExtensives(void) const
{
	return extensives_;
}

vector<Extensive>& hdsim::getAllExtensives(void)
{
	return extensives_;
}

double hdsim::getCellVolume(size_t index) const
{
	return pg_.calcVolume
		(serial_generate(CellEdgesGetter(tess_, static_cast<int>(index))));
}

const CacheData& hdsim::getCacheData(void) const
{
	return cache_data_;
}
