#include <cmath>
#include <algorithm>
#include "hdsim2d.hpp"
#include "hydrodynamics_2d.hpp"
#include "AreaFix.hpp"
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif //RICH_MPI
#include <iostream>
#include "hdf5_diagnostics.hpp"

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
 const CellUpdater& cu
#ifdef RICH_MPI
	,const ProcessorUpdate* proc_update
#endif
	) :
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
  cache_data_(tess, pg)
#ifdef RICH_MPI
	,proc_update_(proc_update)
#endif
{}

hdsim::~hdsim(void) {}

void hdsim::TimeAdvance(void)
{
	vector<Vector2D> point_velocities =
		point_motion_(tess_, cells_, time_);

#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_velocities, true);
	MPI_exchange_data(tess_, cells_, true);
#endif

	vector<Vector2D> edge_velocities =
	  edge_velocity_calculator_(tess_,point_velocities);

	const double dt = tsf_(tess_,
		cells_,
		eos_,
		edge_velocities,
		time_);

	point_velocities=point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_velocities, true);
#endif

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
	if (proc_update_ != 0)
		proc_update_->Update(proctess_, tess_);
	MoveMeshPoints(point_velocities, dt, tess_,proctess_,false);
	// Keep relevant points
	MPI_exchange_data(tess_, extensives_, false);
	MPI_exchange_data(tess_, cells_, false);
#else
	vector<int> HilbertIndeces = MoveMeshPoints(point_velocities, dt, tess_, cycle_ % 25 == 0);
	if (cycle_ % 25 == 0)
	{
		extensives_ = VectorValues(extensives_, HilbertIndeces);
		cells_ = VectorValues(cells_, HilbertIndeces);
		point_velocities = VectorValues(point_velocities, HilbertIndeces);
	}

#endif
	cache_data_.reset();

	cells_ = cu_(tess_, pg_, eos_, extensives_, cells_,
		cache_data_);

	time_ += dt;
	cycle_++;
}

void hdsim::TimeAdvanceClip(void)
{
	vector<Vector2D> point_velocities =
		point_motion_(tess_, cells_, time_);

	vector<Vector2D> edge_velocities =
		edge_velocity_calculator_(tess_, point_velocities);

	const double dt = tsf_(tess_,
		cells_,
		eos_,
		edge_velocities,
		time_);

	point_velocities = point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities);

	edge_velocities =
		edge_velocity_calculator_(tess_, point_velocities);

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

	boost::scoped_ptr<Tessellation> oldtess(tess_.clone());

	MoveMeshPoints(point_velocities, dt, tess_,false);

	extensives_= extensives_+FluxFix2(*oldtess, *oldtess, tess_, point_velocities, dt, cells_, fluxes, edge_velocities, obc_, eos_);

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

#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_velocities, true);
	MPI_exchange_data(tess_, cells_, true);
#endif


	vector<Vector2D> edge_velocities =
	  edge_velocity_calculator_(tess_,point_velocities);

	const double dt = tsf_(tess_, cells_, eos_, edge_velocities, time_);

	point_velocities=point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities);

#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_velocities, true);
#endif
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
	
	eu_(mid_fluxes, pg_, tess_, dt, cache_data_, cells_, mid_extensives);


#ifdef RICH_MPI
	if (proc_update_ != 0)
		proc_update_->Update(proctess_, tess_);
	vector<int> HilbertIndeces = MoveMeshPoints(point_velocities, dt, tess_, proctess_,cycle_%25==0);
	if (cycle_ % 25 == 0)
	{
		mid_extensives = VectorValues(mid_extensives, HilbertIndeces);
		extensives_ = VectorValues(extensives_, HilbertIndeces);
		cells_ = VectorValues(cells_, HilbertIndeces);
		point_velocities = VectorValues(point_velocities, HilbertIndeces);
	}
	// Keep relevant points
	MPI_exchange_data(tess_, mid_extensives, false);
	MPI_exchange_data(tess_, extensives_, false);
	MPI_exchange_data(tess_, cells_, false);
	MPI_exchange_data(tess_, point_velocities, false);
	MPI_exchange_data(tess_, point_velocities, true);
#else
	vector<int> HilbertIndeces = MoveMeshPoints(point_velocities, dt, tess_,cycle_%25==0);
	if (cycle_ % 25 == 0)
	{
		mid_extensives = VectorValues(mid_extensives, HilbertIndeces);
		extensives_ = VectorValues(extensives_, HilbertIndeces);
		cells_ = VectorValues(cells_, HilbertIndeces);
		point_velocities = VectorValues(point_velocities, HilbertIndeces);
	}
#endif
	cache_data_.reset();

	vector<ComputationalCell> mid_cells = cu_(tess_, pg_, eos_, mid_extensives, cells_, cache_data_);

#ifdef RICH_MPI
	MPI_exchange_data(tess_, mid_cells,true);
#endif
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
	   mid_extensives);

	eu_(fluxes, pg_, tess_, dt, cache_data_, cells_, mid_extensives);


	extensives_ = average_extensive(mid_extensives, extensives_);

	cells_ = cu_(tess_, pg_, eos_, extensives_, cells_, cache_data_);

	time_ += dt;
	++cycle_;
}

void hdsim::TimeAdvance2MidPointClip(void)
{
	vector<Vector2D> point_velocities = point_motion_(tess_, cells_, time_);

	vector<Vector2D> edge_velocities =edge_velocity_calculator_(tess_, point_velocities);

	const double dt = tsf_(tess_, cells_, eos_, edge_velocities, time_);

	point_velocities = point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities);

	edge_velocities = edge_velocity_calculator_(tess_, point_velocities);

	vector<Extensive> fluxes = fc_(tess_,edge_velocities,cells_,extensives_,cache_data_,eos_,time_,0.5*dt);

	vector<Extensive> mid_extensives = extensives_;

	boost::scoped_ptr<Tessellation> oldtess(tess_.clone());

	vector<Vector2D> old_points = tess_.GetMeshPoints();
	old_points.resize(static_cast<size_t>(tess_.GetPointNo()));
	MoveMeshPoints(point_velocities, 0.5*dt, tess_,false);

	CacheData data_temp(*oldtess, pg_);

	mid_extensives = mid_extensives + FluxFix2(*oldtess, *oldtess, tess_, point_velocities, 0.5*dt, cells_, fluxes,
		edge_velocities, obc_, eos_);

	eu_(fluxes, pg_, *oldtess, 0.5*dt, data_temp, cells_, mid_extensives);

	ExternalForceContribution(*oldtess,pg_, data_temp,cells_,fluxes,point_velocities,source_,time_,0.5*dt,mid_extensives);

	time_ += 0.5*dt;
	cache_data_.reset();

	vector<ComputationalCell> mid_cells = cu_(tess_, pg_, eos_, mid_extensives, cells_, cache_data_);

	edge_velocities = edge_velocity_calculator_(tess_, point_velocities);

	fluxes = fc_(tess_,edge_velocities,mid_cells,mid_extensives,cache_data_,eos_,time_,dt);

	boost::scoped_ptr<Tessellation> midtess(tess_.clone());

	MoveMeshPoints(point_velocities, dt, tess_,false, old_points);

	CacheData cachetemp2(*midtess, pg_);

	extensives_ = extensives_ + FluxFix2(*oldtess, *midtess, tess_, point_velocities, dt, mid_cells, fluxes,
		edge_velocities, obc_, eos_);

	eu_(fluxes, pg_, *midtess, dt, cachetemp2, cells_, extensives_);

	ExternalForceContribution(*midtess,pg_, cachetemp2,mid_cells,fluxes,point_velocities,source_,time_,dt,extensives_);

	cache_data_.reset();
	cells_ = cu_(tess_, pg_, eos_, extensives_, cells_, cache_data_);

	if (cycle_ % 25 == 0)
	{
		vector<int> HilbertIndeces = MoveMeshPoints(point_velocities, dt, tess_, cycle_ % 25 == 0,old_points);
		extensives_ = VectorValues(extensives_, HilbertIndeces);
		cells_ = VectorValues(cells_, HilbertIndeces);
		cache_data_.reset();
	}

	time_ += 0.5*dt;
	++cycle_;
}

void hdsim::TimeAdvance2MidPoint(void)
{
	vector<Vector2D> point_velocities = point_motion_(tess_, cells_, time_);

	vector<Vector2D> edge_velocities = edge_velocity_calculator_(tess_, point_velocities);

	const double dt = tsf_(tess_, cells_, eos_, edge_velocities, time_);

	point_velocities = point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities);

	edge_velocities = edge_velocity_calculator_(tess_, point_velocities);

	vector<Extensive> fluxes = fc_(tess_, edge_velocities, cells_, extensives_, cache_data_, eos_, time_, 0.5*dt);

	vector<Extensive> mid_extensives = extensives_;

	eu_(fluxes, pg_, tess_, 0.5*dt, cache_data_, cells_, mid_extensives);

	ExternalForceContribution(tess_, pg_, cache_data_, cells_, fluxes, point_velocities, source_, time_, 0.5*dt, mid_extensives);

	vector<Vector2D> old_points = tess_.GetMeshPoints();
	old_points.resize(static_cast<size_t>(tess_.GetPointNo()));

	boost::scoped_ptr<Tessellation> oldtess(tess_.clone());

	vector<int> HilbertIndeces = MoveMeshPoints(point_velocities, 0.5*dt, tess_, cycle_ % 25 == 0);
	if (cycle_ % 25 == 0)
	{
		mid_extensives = VectorValues(mid_extensives, HilbertIndeces);
		extensives_ = VectorValues(extensives_, HilbertIndeces);
		cells_ = VectorValues(cells_, HilbertIndeces);
		point_velocities = VectorValues(point_velocities, HilbertIndeces);
		old_points = VectorValues(old_points, HilbertIndeces);
	}


	time_ += 0.5*dt;

	cache_data_.reset();

	vector<ComputationalCell> mid_cells = cu_(tess_, pg_, eos_, mid_extensives, cells_, cache_data_);

	edge_velocities = edge_velocity_calculator_(tess_, point_velocities);

	fluxes = fc_(tess_, edge_velocities, mid_cells, mid_extensives, cache_data_, eos_, time_, dt);

	eu_(fluxes, pg_, tess_, dt, cache_data_, cells_, extensives_);

	ExternalForceContribution(tess_, pg_, cache_data_, mid_cells, fluxes, point_velocities, source_, time_, dt, extensives_);

	boost::scoped_ptr<Tessellation> midtess(tess_.clone());

	MoveMeshPoints(point_velocities, dt, tess_, false,old_points);
	cache_data_.reset();

	cells_ = cu_(tess_, pg_, eos_, extensives_, cells_, cache_data_);

	time_ += 0.5*dt;
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
		if (extensives_[0].tracers.size() != cell.tracers.size())
		{
			for (boost::container::flat_map<std::string, double>::const_iterator it =
				cell.tracers.begin();
				it != cell.tracers.end(); ++it)
				extensives_[i].tracers[it->first] = (it->second)*mass;
		}
		else
		{
			boost::container::flat_map<string, double>::const_iterator it2 = cell.tracers.begin();
			for (boost::container::flat_map<string, double>::iterator it = extensives_[i].tracers.begin(); it !=
				extensives_[i].tracers.end(); ++it, ++it2)
				it->second = it2->second*mass;
		}
	}
}

void hdsim::setCycle(int cycle)
{
	cycle_ = cycle;
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

#ifdef RICH_MPI
const Tessellation & hdsim::GetProcTessellation(void)const
{
	return proctess_;
}
#endif// RICH_MPI
