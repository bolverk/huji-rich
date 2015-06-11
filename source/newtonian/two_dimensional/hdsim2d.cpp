#include <cmath>
#include <algorithm>
#include "hdsim2d.hpp"
#include "hydrodynamics_2d.hpp"

using namespace std;

namespace {
  class CellEdgesGetter: public LazyList<Edge>
  {
  public:

    CellEdgesGetter(const Tessellation& tess, int n):
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

namespace {

  vector<Extensive> init_extensives(const Tessellation& tess,
				    const PhysicalGeometry& pg,
				    const vector<ComputationalCell>& cells,
				    const EquationOfState& eos)
  {
    vector<Extensive> res(cells.size());
    for(size_t i=0;i<cells.size();++i){
      const ComputationalCell& cell = cells[i];
      const double volume = 
	pg.calcVolume
	(serial_generate(CellEdgesGetter(tess,static_cast<int>(i))));
      const double mass = volume*cell.density;
      res[i].mass = mass;
      res[i].energy = eos.dp2e(cell.density, cell.pressure,cell.tracers)*mass + 
	0.5*mass*ScalarProd(cell.velocity, cell.velocity);
      res[i].momentum = mass*cells[i].velocity;
      for(std::map<std::string,double>::const_iterator it=
	    cells[i].tracers.begin();
	  it!=cells[i].tracers.end();++it)
	res[i].tracers[it->first] = (it->second)*mass;
    }
    return res;
  }
}

hdsim::hdsim
(Tessellation& tess,
 const OuterBoundary& obc,
 const PhysicalGeometry& pg,
 const vector<ComputationalCell>& cells,
 const EquationOfState& eos,
 const PointMotion& point_motion,
 const SourceTerm& source,
 const TimeStepFunction& tsf,
 const FluxCalculator& fc,
 const ExtensiveUpdater& eu,
 const CellUpdater& cu):
  tess_(tess),
  obc_(obc),
  eos_(eos),
  cells_(cells),
  extensives_(init_extensives(tess,
			      pg,
			      cells,
			      eos)),
  point_motion_(point_motion),
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

namespace
{
  /*
    void update_extensives(const vector<Extensive>& fluxes,
    const PhysicalGeometry& pg,
    const Tessellation& tess,
    const double dt,
    vector<Extensive>& extensives)
    {
    const vector<Edge>& edge_list = tess.getAllEdges();
    for(size_t i=0;i<edge_list.size();++i){
    const Edge edge = edge_list[i];
    const Extensive delta = dt*pg.calcArea(edge)*fluxes[i];
    if(edge.neighbors.first>=0 && 
    edge.neighbors.first < static_cast<int>(extensives.size()))
    extensives[static_cast<size_t>(edge.neighbors.first)] -= delta;
    if(edge.neighbors.second>=0 &&
    edge.neighbors.second < static_cast<int>(extensives.size()))
    extensives[static_cast<size_t>(edge.neighbors.second)] += delta;
    }
    }
  */
}

void hdsim::TimeAdvance(void)
{
  const vector<Vector2D> point_velocities = 
    point_motion_(tess_,cells_,time_);

  const double dt = tsf_(tess_,
			 cells_,
			 eos_,
			 point_velocities,
			 time_);
			 
  const vector<Extensive> fluxes = fc_(tess_,
				       point_velocities,
				       cells_,
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

  MoveMeshPoints(point_velocities, dt, tess_);
  cache_data_.reset();

  cells_ = cu_(tess_, pg_, eos_, extensives_, cells_,
	       cache_data_);

  time_ += dt;
  cycle_++;
}

namespace {

  template<class T> class AverageCalculator: public LazyList<T>
  {
  public:

    AverageCalculator(const vector<T>& ll1,
		      const vector<T>& ll2):
      ll1_(ll1), ll2_(ll2)
    {
      assert(ll1.size()==ll2.size());
    }

    size_t getLength(void) const
    {
      return ll1_.size();
    }

    T operator()(size_t i) const
    {
      return 0.5*(ll1_[i]+ll2_[i]);
    }

  private:
    const vector<T>& ll1_;
    const vector<T>& ll2_;
  };

  template<class T> vector<T> average(const vector<T>& v1,
				      const vector<T>& v2)
  {
    return serial_generate(AverageCalculator<T>(v1,v2));
  }
}

const PhysicalGeometry& hdsim::getPhysicalGeometry(void) const
{
  return pg_;
}

const Tessellation& hdsim::getTessellation(void) const
{
  return tess_;
}

#ifdef RICH_MPI
Tessellation const& hdsim::GetProcTessellation(void) const
{
  return _proctess;
}
#endif

void hdsim::addTracer(const string& name,
		      const SpatialDistribution& tp)
{
  for(size_t i=0;i<cells_.size();++i)
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

namespace
{
#ifdef RICH_MPI
  void CreateGetPrimitiveList(vector<int> const& ToRemove,vector<vector<int> >
			      const& Nghost,int /*nremoved*/,vector<vector<int> > &MPI_AMR_Send)
  {
    int nprocs=static_cast<int>(Nghost.size());
    MPI_AMR_Send.resize(static_cast<size_t>(nprocs));
    vector<vector<int> > SortedNghost(static_cast<size_t>(nprocs)),SortIndex(static_cast<size_t>(nprocs));
    // sort Nghost
    for(size_t i=0;i<static_cast<size_t>(nprocs);++i)
      {
	SortedNghost[i]=Nghost[i];
	sort_index(Nghost[i],SortIndex[i]);
	sort(SortedNghost[i].begin(),SortedNghost[i].end());
      }
    for(int i=0;i<static_cast<int>(ToRemove.size());++i)
      {
	for(int j=0;j<nprocs;++j)
	  {
	    if(binary_search(SortedNghost[static_cast<size_t>(j)].begin(),SortedNghost[static_cast<size_t>(j)].end(),
			     ToRemove[static_cast<size_t>(i)]))
	      {
		int index2=static_cast<int>(lower_bound(SortedNghost[static_cast<size_t>(j)].begin(),SortedNghost[static_cast<size_t>(j)].end(),
							ToRemove[static_cast<size_t>(i)])-SortedNghost[static_cast<size_t>(j)].begin());
		MPI_AMR_Send[static_cast<size_t>(j)].push_back(SortIndex[static_cast<size_t>(j)][static_cast<size_t>(index2)]);
	      }
	  }
      }
  }
#endif // RICH_MPI
}

namespace
{
#ifdef RICH_MPI
  void RemoveNGhostAMR(vector<vector<int> > &nghost,vector<int> const& sentprocs,
		       vector<vector<int> > &toremove)
  {
    int nlist=static_cast<int>(sentprocs.size());
    int rank=get_mpi_rank();
    int ws=get_mpi_size();
    vector<int> procorder=GetProcOrder(rank,ws);
    vector<vector<int> > recv(static_cast<size_t>(nlist));
    int temp;
    MPI_Status status;
    for(int i=0;i<static_cast<int>(procorder.size());++i)
      {
	int index=static_cast<int>(Find(sentprocs.begin(),sentprocs.end(),procorder[static_cast<size_t>(i)])
				   -sentprocs.begin());
	if(index<nlist)
	  {
	    if(rank<procorder[static_cast<size_t>(i)])
	      {
		if(toremove[static_cast<size_t>(index)].empty())
		  MPI_Send(&temp,1,MPI_INT,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD);
		else
		  MPI_Send(&toremove[static_cast<size_t>(index)][0],static_cast<int>(toremove[static_cast<size_t>(index)].size()),
			   MPI_INT,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD);
		MPI_Probe(procorder[static_cast<size_t>(i)],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		if(status.MPI_TAG==1)
		  MPI_Recv(&temp,1,MPI_INT,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		else
		  {
		    int count;
		    MPI_Get_count(&status,MPI_INT,&count);
		    recv[static_cast<size_t>(index)].resize(static_cast<size_t>(count));
		    MPI_Recv(&recv[static_cast<size_t>(index)][0],count,MPI_INT,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		  }
	      }
	    else
	      {
		MPI_Probe(procorder[static_cast<size_t>(i)],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		if(status.MPI_TAG==1)
		  MPI_Recv(&temp,1,MPI_INT,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		else
		  {
		    int count;
		    MPI_Get_count(&status,MPI_INT,&count);
		    recv[static_cast<size_t>(index)].resize(static_cast<size_t>(count));
		    MPI_Recv(&recv[static_cast<size_t>(index)][0],count,MPI_INT,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		  }
		if(toremove[static_cast<size_t>(index)].empty())
		  MPI_Send(&temp,1,MPI_INT,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD);
		else
		  MPI_Send(&toremove[static_cast<size_t>(index)][0],static_cast<int>(toremove[static_cast<size_t>(index)].size()),
			   MPI_INT,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD);
	      }
	  }
      }
    for(int i=0;i<nlist;++i)
      if(!recv[static_cast<size_t>(i)].empty())
	RemoveVector(nghost[static_cast<size_t>(i)],recv[static_cast<size_t>(i)]);
  }
#endif
}

const vector<ComputationalCell>& hdsim::getAllCells(void) const
{
  return cells_;
}

vector<ComputationalCell>& hdsim::getAllCells(void)
{
  return cells_;
}

void hdsim::recalculateExtensives(void)
{
  for(size_t i=0;i<extensives_.size();++i){
    const ComputationalCell& cell = cells_[i];
    const double volume = cache_data_.volumes[i];
    const double mass = volume*cell.density;
    extensives_[i].mass = mass;
    extensives_[i].energy = eos_.dp2e(cell.density, cell.pressure)*mass +
      0.5*mass*ScalarProd(cell.velocity, cell.velocity);
    extensives_[i].momentum = mass*cell.velocity;
    for(std::map<std::string,double>::const_iterator it=
	  cell.tracers.begin();
	it!=cell.tracers.end();++it)
      extensives_[i].tracers[it->first] = (it->second)*mass;
  }
}

/*
  void hdsim::HilbertArrange(int innernum)
  {
  vector<Vector2D> cor=_tessellation.GetMeshPoints();
  vector<int> order=HilbertOrder(cor,_tessellation.GetPointNo(),innernum);
  ReArrangeVector(cor,order);
  if(cor.size()>order.size())
  cor.erase(cor.begin()+static_cast<int>(order.size()),cor.end());
  ReArrangeVector(_cells,order);
  if(tracer_flag_)
  ReArrangeVector(tracer_,order);
  _tessellation.Update(cor);
  }
*/

#ifdef RICH_MPI
void hdsim::SetProcessorMovement(ProcessorUpdate *procupdate)
{
  procupdate_=procupdate;
}
#endif

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

double hdsim::getCellVolume(size_t index) const
{
  return pg_.calcVolume
    (serial_generate(CellEdgesGetter(tess_,static_cast<int>(index))));
}
