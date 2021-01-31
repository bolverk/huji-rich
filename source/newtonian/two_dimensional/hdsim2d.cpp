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
#include "../../misc/serial_generate.hpp"
#include <numeric>


using namespace std;

namespace
{

  Extensive cell2extensive
  (const ComputationalCell& cell,
   const double volume,
   const EquationOfState& eos,
   const TracerStickerNames& tnames)
  {
    const double mass = volume*cell.density;
    Extensive res;
    res.mass = mass;
    res.energy = eos.dp2e
      (cell.density,
       cell.pressure,
       cell.tracers,
       tnames.tracer_names)*mass+
      0.5*mass*ScalarProd(cell.velocity,
			  cell.velocity);
    res.momentum = mass*cell.velocity;
    res.tracers = serial_generate<double, double>
      (cell.tracers,
       [&mass](double t){return mass*t;});
    return res;
  }

  Extensive relativistic_cell2extensive
  (const ComputationalCell& cell,
   const double volume,
   const EquationOfState& eos,
   const TracerStickerNames& tnames)
  {
    const double gamma = 1 / std::sqrt(1 - ScalarProd(cell.velocity, cell.velocity));
    const double mass = volume*cell.density*gamma;
    const double enthalpy = eos.dp2e(cell.density,
				     cell.pressure,
				     cell.tracers,
				     tnames.tracer_names);
    Extensive res;
    res.mass = mass;
    if (fastabs(cell.velocity) < 1e-5)
      res.energy = (gamma*enthalpy + 0.5*ScalarProd(cell.velocity, cell.velocity))* mass - cell.pressure*volume;
    else
      res.energy = (gamma*enthalpy + (gamma - 1))* mass - cell.pressure*volume;
    res.momentum = mass * (enthalpy+1)*gamma*cell.velocity;
    res.tracers = serial_generate<double, double>
      (cell.tracers,
       [&mass](double t){return mass*t;});
    return res;
  }

  namespace {
    vector<size_t> create_range(size_t n)
    {
      vector<size_t> res(n);
      std::iota(res.begin(),
		res.end(),
		0);
      return res;
    }
  }
  
  vector<Extensive> init_extensives(const Tessellation& tess,
				    const PhysicalGeometry& pg,
				    const vector<ComputationalCell>& cells,
				    const EquationOfState& eos,
				    TracerStickerNames const& tracernames,
				    bool relativistic)
  {
    return serial_generate<size_t, Extensive>
      (create_range(static_cast<size_t>(tess.GetPointNo())),
       [&](size_t i){
	 const ComputationalCell& cell = cells[i];
	 const double volume =
	   pg.calcVolume
	   (serial_generate<int, Edge>
	    (tess.GetCellEdges(static_cast<int>(i)),
	     [&tess](int j){return tess.GetEdge(j);}));
	 auto func = (relativistic ?
		      &relativistic_cell2extensive :
		      &cell2extensive);
	 return (*func)
	   (cell,
	    volume,
	    eos,
	    tracernames);
       });
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
 const CellUpdater& cu,
 const TracerStickerNames& tracer_sticker_names,
 bool relativistic
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
    eos,
    tracer_sticker_names,relativistic)),
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
  tracer_sticker_names_(tracer_sticker_names),
  cache_data_(tess, pg)
#ifdef RICH_MPI
  , proc_update_(proc_update)
#endif
{
  // sort tracers and stickers
  const vector<size_t> tindex =
    sort_index(tracer_sticker_names_.tracer_names);
  const vector<size_t> sindex =
    sort_index(tracer_sticker_names_.sticker_names);
  tracer_sticker_names_.tracer_names =
    VectorValues(tracer_sticker_names_.tracer_names,
		 tindex);
  tracer_sticker_names_.sticker_names =
    VectorValues(tracer_sticker_names_.sticker_names,
		 sindex);
  for_each(cells_.begin(),
	   cells_.end(),
	   [&](ComputationalCell& c){
	     c.tracers = serial_generate<size_t, double>
	       (tindex,
		[&](size_t j){return c.tracers.at(j);});
	   });
  for(pair<vector<Extensive>::iterator, vector<ComputationalCell>::iterator > multitr{extensives_.begin(), cells_.begin()};
      multitr.first!=extensives_.end();
      ++multitr.first,++multitr.second)
    (*multitr.first).tracers = serial_generate<double, double>
      ((*multitr.second).tracers,
       [&](double t){return t*(*multitr.first).mass;});
  for_each(cells_.begin(),
	   cells_.end(),
	   [&](ComputationalCell& c){
	     c.stickers = serial_generate<size_t, bool>
	       (sindex,
		[&](size_t j){return c.stickers.at(j);});});
#ifdef RICH_MPI
  MPI_exchange_data(tess_, cells_, true);
#endif
}

hdsim::~hdsim(void) {}

TracerStickerNames const& hdsim::GetTracerStickerNames(void)const
{
  return tracer_sticker_names_;
}

void hdsim::TimeAdvance(void)
{
  vector<Vector2D> point_velocities =
    point_motion_(tess_, cells_, time_,tracer_sticker_names_);

#ifdef RICH_MPI
  MPI_exchange_data(tess_, point_velocities, true);
#endif

  vector<Vector2D> edge_velocities =
    edge_velocity_calculator_(tess_, point_velocities);

  const double dt = tsf_(tess_,
			 cells_,
			 eos_,
			 edge_velocities,
			 time_,tracer_sticker_names_);

  point_velocities = point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities,tracer_sticker_names_);
#ifdef RICH_MPI
  MPI_exchange_data(tess_, point_velocities, true);
#endif

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
     dt,tracer_sticker_names_);


  //  update_extensives(fluxes,
  eu_(fluxes,
      pg_,
      tess_,
      dt,
      cache_data_,
      cells_,
      extensives_,
      time_,tracer_sticker_names_);

  ExternalForceContribution(tess_,
			    pg_,
			    cache_data_,
			    cells_,
			    fluxes,
			    point_velocities,
			    source_,
			    time_,
			    dt,
			    extensives_,tracer_sticker_names_);

#ifdef RICH_MPI
  if (proc_update_ != 0)
    proc_update_->Update(proctess_, tess_);
  vector<int> HilbertIndeces = MoveMeshPoints(point_velocities, dt, tess_, proctess_, cycle_ % 10 == 0);
  if (cycle_ % 10 == 0)
    {
      extensives_ = VectorValues(extensives_, HilbertIndeces);
      cells_ = VectorValues(cells_, HilbertIndeces);
      point_velocities = VectorValues(point_velocities, HilbertIndeces);
    }
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

  cells_ = cu_(tess_, pg_, eos_, extensives_, cells_,cache_data_,tracer_sticker_names_, time_);

#ifdef RICH_MPI
  MPI_exchange_data(tess_, cells_, true);
#endif

  time_ += dt;
  cycle_++;
}

void hdsim::TimeAdvanceClip(void)
{
  vector<Vector2D> point_velocities =
    point_motion_(tess_, cells_, time_,tracer_sticker_names_);

  vector<Vector2D> edge_velocities =
    edge_velocity_calculator_(tess_, point_velocities);

  const double dt = tsf_(tess_,
			 cells_,
			 eos_,
			 edge_velocities,
			 time_,tracer_sticker_names_);

  point_velocities = point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities,tracer_sticker_names_);

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
     dt,tracer_sticker_names_);


  //  update_extensives(fluxes,
  eu_(fluxes,
      pg_,
      tess_,
      dt,
      cache_data_,
      cells_,
      extensives_, time_,tracer_sticker_names_);

  ExternalForceContribution(tess_,
			    pg_,
			    cache_data_,
			    cells_,
			    fluxes,
			    point_velocities,
			    source_,
			    time_,
			    dt,
			    extensives_,tracer_sticker_names_);
	
  boost::scoped_ptr<Tessellation> oldtess(tess_.clone());

  MoveMeshPoints(point_velocities, dt, tess_, false);

  extensives_ = extensives_ + FluxFix2(*oldtess, *oldtess, tess_, point_velocities, dt, cells_, fluxes, edge_velocities, obc_, eos_);

  cache_data_.reset();

  cells_ = cu_(tess_, pg_, eos_, extensives_, cells_,cache_data_,tracer_sticker_names_, time_);

  time_ += dt;
  cycle_++;
}

void hdsim::TimeAdvance2Heun(void)
{
#ifdef RICH_DEBUG_PRINT
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "Here 0" << std::endl;
#endif
  vector<Vector2D> point_velocities = point_motion_(tess_, cells_, time_,tracer_sticker_names_);

#ifdef RICH_MPI
  MPI_exchange_data(tess_, point_velocities, true);
#endif

  vector<Vector2D> edge_velocities =
    edge_velocity_calculator_(tess_, point_velocities);
#ifdef RICH_DEBUG_PRINT
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "Here 1" << std::endl;
#endif

  double dt = tsf_(tess_, cells_, eos_, edge_velocities, time_,tracer_sticker_names_);

#ifdef RICH_DEBUG_PRINT
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "Here 2" << std::endl;
#endif

  point_velocities = point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities,tracer_sticker_names_);

#ifdef RICH_MPI
  MPI_exchange_data(tess_, point_velocities, true);
#endif
  edge_velocities =
    edge_velocity_calculator_(tess_, point_velocities);
	
  dt = tsf_(tess_, cells_, eos_, edge_velocities, time_, tracer_sticker_names_);
#ifdef RICH_DEBUG_PRINT
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "Here 3" << std::endl;
#endif

  const vector<Extensive> mid_fluxes =
    fc_
    (tess_,
     edge_velocities,
     cells_,
     extensives_,
     cache_data_,
     eos_,
     time_,
     dt,tracer_sticker_names_);
#ifdef RICH_DEBUG_PRINT
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "Here 4" << std::endl;
#endif

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
     mid_extensives,tracer_sticker_names_);
#ifdef RICH_DEBUG_PRINT
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "Here 5" << std::endl;
#endif

  eu_(mid_fluxes, pg_, tess_, dt, cache_data_, cells_, mid_extensives, time_,tracer_sticker_names_);

#ifdef RICH_DEBUG_PRINT
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "Here 6" << std::endl;
#endif

#ifdef RICH_MPI
  if (proc_update_ != 0)
    proc_update_->Update(proctess_, tess_);
  vector<int> HilbertIndeces = MoveMeshPoints(point_velocities, dt, tess_, proctess_, cycle_ % 10 == 0);
  if (cycle_ % 10 == 0)
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
  vector<int> HilbertIndeces = MoveMeshPoints(point_velocities, dt, tess_, cycle_ % 25 == 0);
  if (cycle_ % 25 == 0)
    {
      mid_extensives = VectorValues(mid_extensives, HilbertIndeces);
      extensives_ = VectorValues(extensives_, HilbertIndeces);
      cells_ = VectorValues(cells_, HilbertIndeces);
      point_velocities = VectorValues(point_velocities, HilbertIndeces);
    }
#endif
#ifdef RICH_DEBUG_PRINT
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "Here 7" << std::endl;
#endif

  cache_data_.reset();

  vector<ComputationalCell> mid_cells = cu_(tess_, pg_, eos_, mid_extensives, cells_, cache_data_,
					    tracer_sticker_names_, time_ + dt);
#ifdef RICH_DEBUG_PRINT
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "Here 8" << std::endl;
#endif

#ifdef RICH_MPI
  MPI_exchange_data(tess_, mid_cells, true);
  MPI_exchange_data(tess_, cells_, true);
#endif
  edge_velocities =
    edge_velocity_calculator_(tess_, point_velocities);

  const vector<Extensive> fluxes =
    fc_
    (tess_,
     edge_velocities,
     mid_cells,
     mid_extensives,
     cache_data_,
     eos_,
     time_,
     dt,tracer_sticker_names_);

#ifdef RICH_DEBUG_PRINT
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "Here 9" << std::endl;
#endif

  ExternalForceContribution
    (tess_,
     pg_,
     cache_data_,
     mid_cells,
     fluxes,
     point_velocities,
     source_,
     time_ + dt,
     dt,
     mid_extensives,tracer_sticker_names_);

  eu_(fluxes, pg_, tess_, dt, cache_data_, cells_, mid_extensives, time_ + dt,tracer_sticker_names_);

  transform(extensives_.begin(), extensives_.end(),
	    mid_extensives.begin(), extensives_.begin(),
	    [](const Extensive& arg1, const Extensive& arg2){
	      return 0.5*(arg1+arg2);
	    });
	    

  cells_ = cu_(tess_, pg_, eos_, extensives_, cells_, cache_data_,tracer_sticker_names_, time_ + dt);

#ifdef RICH_DEBUG_PRINT
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "Here 11" << std::endl;
#endif

#ifdef RICH_MPI
  MPI_exchange_data(tess_, cells_, true);
#endif

  time_ += dt;
  ++cycle_;
}

void hdsim::TimeAdvance2MidPointClip(void)
{
  vector<Vector2D> point_velocities = point_motion_(tess_, cells_, time_,tracer_sticker_names_);

  vector<Vector2D> edge_velocities = edge_velocity_calculator_(tess_, point_velocities);

  const double dt = tsf_(tess_, cells_, eos_, edge_velocities, time_,tracer_sticker_names_);

  point_velocities = point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities,tracer_sticker_names_);

  edge_velocities = edge_velocity_calculator_(tess_, point_velocities);

  vector<Extensive> fluxes = fc_(tess_, edge_velocities, cells_, extensives_, cache_data_, eos_, time_, 0.5*dt,
				 tracer_sticker_names_);

  vector<Extensive> mid_extensives = extensives_;

  boost::scoped_ptr<Tessellation> oldtess(tess_.clone());

  vector<Vector2D> old_points = tess_.GetMeshPoints();
  old_points.resize(static_cast<size_t>(tess_.GetPointNo()));
  MoveMeshPoints(point_velocities, 0.5*dt, tess_, false);

  CacheData data_temp(*oldtess, pg_);

  mid_extensives = mid_extensives + FluxFix2(*oldtess, *oldtess, tess_, point_velocities, 0.5*dt, cells_, fluxes,
					     edge_velocities, obc_, eos_);

  eu_(fluxes, pg_, *oldtess, 0.5*dt, data_temp, cells_, mid_extensives, time_,tracer_sticker_names_);

  ExternalForceContribution(*oldtess, pg_, data_temp, cells_, fluxes, point_velocities, source_, time_, 0.5*dt,
			    mid_extensives,tracer_sticker_names_);

  time_ += 0.5*dt;
  cache_data_.reset();

  vector<ComputationalCell> mid_cells = cu_(tess_, pg_, eos_, mid_extensives, cells_, cache_data_,tracer_sticker_names_, time_);

  edge_velocities = edge_velocity_calculator_(tess_, point_velocities);

  fluxes = fc_(tess_, edge_velocities, mid_cells, mid_extensives, cache_data_, eos_, time_, dt,tracer_sticker_names_);

  boost::scoped_ptr<Tessellation> midtess(tess_.clone());

  MoveMeshPoints(point_velocities, dt, tess_, false, old_points);

  CacheData cachetemp2(*midtess, pg_);

  extensives_ = extensives_ + FluxFix2(*oldtess, *midtess, tess_, point_velocities, dt, mid_cells, fluxes,
				       edge_velocities, obc_, eos_);

  eu_(fluxes, pg_, *midtess, dt, cachetemp2, cells_, extensives_, time_,tracer_sticker_names_);

  ExternalForceContribution(*midtess, pg_, cachetemp2, mid_cells, fluxes, point_velocities, source_, time_, dt, 
			    extensives_,tracer_sticker_names_);

  cache_data_.reset();
  cells_ = cu_(tess_, pg_, eos_, extensives_, cells_, cache_data_,tracer_sticker_names_, time_ + 0.5 * dt);

  if (cycle_ % 25 == 0)
    {
      vector<int> HilbertIndeces = MoveMeshPoints(point_velocities, dt, tess_, cycle_ % 25 == 0, old_points);
      extensives_ = VectorValues(extensives_, HilbertIndeces);
      cells_ = VectorValues(cells_, HilbertIndeces);
      cache_data_.reset();
    }

  time_ += 0.5*dt;
  ++cycle_;
}

void hdsim::TimeAdvance2MidPoint(void)
{
  vector<Vector2D> point_velocities = point_motion_(tess_, cells_, time_,tracer_sticker_names_);

  vector<Vector2D> edge_velocities = edge_velocity_calculator_(tess_, point_velocities);

  const double dt = tsf_(tess_, cells_, eos_, edge_velocities, time_,tracer_sticker_names_);

  point_velocities = point_motion_.ApplyFix(tess_, cells_, time_, dt, point_velocities,tracer_sticker_names_);

  edge_velocities = edge_velocity_calculator_(tess_, point_velocities);

  vector<Extensive> fluxes = fc_(tess_, edge_velocities, cells_, extensives_, cache_data_, eos_, time_, 0.5*dt,
				 tracer_sticker_names_);

  vector<Extensive> mid_extensives = extensives_;

  eu_(fluxes, pg_, tess_, 0.5*dt, cache_data_, cells_, mid_extensives, time_,tracer_sticker_names_);

  ExternalForceContribution(tess_, pg_, cache_data_, cells_, fluxes, point_velocities, source_, time_, 0.5*dt, 
			    mid_extensives,tracer_sticker_names_);

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

  vector<ComputationalCell> mid_cells = cu_(tess_, pg_, eos_, mid_extensives, cells_, cache_data_,
					    tracer_sticker_names_, time_);

  edge_velocities = edge_velocity_calculator_(tess_, point_velocities);

  fluxes = fc_(tess_, edge_velocities, mid_cells, mid_extensives, cache_data_, eos_, time_, dt,
	       tracer_sticker_names_);

  eu_(fluxes, pg_, tess_, dt, cache_data_, cells_, extensives_, time_,tracer_sticker_names_);

  ExternalForceContribution(tess_, pg_, cache_data_, mid_cells, fluxes, point_velocities, source_, time_, dt, 
			    extensives_,tracer_sticker_names_);

  boost::scoped_ptr<Tessellation> midtess(tess_.clone());

  MoveMeshPoints(point_velocities, dt, tess_, false, old_points);
  cache_data_.reset();

  cells_ = cu_(tess_, pg_, eos_, extensives_, cells_, cache_data_,tracer_sticker_names_, time_ + 0.5 * dt);

  time_ += 0.5*dt;
  ++cycle_;
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
	       cache_data_,tracer_sticker_names_, time_);
#ifdef RICH_MPI
  MPI_exchange_data(tess_, cells_, true);
#endif

}

void hdsim::recalculateExtensives(void)
{
  const vector<double> volumes =
    serial_generate<size_t, double>
    (create_range(static_cast<size_t>(tess_.GetPointNo())),
     [&](size_t i){return cache_data_.volumes[i];});
  transform(cells_.begin(),
	    cells_.end(),
	    volumes.begin(),
	    extensives_.begin(),
	    [&](const ComputationalCell& c,
		const double& v){
	      return cell2extensive(c,
				    v,
				    eos_,
				    tracer_sticker_names_);
	    });
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
    (serial_generate<int, Edge>
     (tess_.GetCellEdges(static_cast<int>(index)),
      [&](int j){return tess_.GetEdge(j);}));
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
