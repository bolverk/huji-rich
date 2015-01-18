#include <cmath>
#include <algorithm>
#include "hdsim2d.hpp"
#include "hydrodynamics_2d.hpp"

using namespace std;

ColdflowParams::ColdflowParams(const double as_i,
			       const double bs_i):
  as(as_i), bs(bs_i) {}

void hdsim::SetData(vector<Primitive> const& cells,vector<Vector2D> const& points,
	double time,vector<vector<double> > const& tracers)
{
	_cells=cells;
	_tessellation.Update(points);
	tracer_=tracers;
	_time=time;
}

hdsim::hdsim
	(vector<Vector2D> const& points,
	Tessellation& tessellation,
#ifdef RICH_MPI
	Tessellation& proctess,
#endif
	SpatialReconstruction& interpolation,
	SpatialDistribution const& density,
	SpatialDistribution const& pressure,
	SpatialDistribution const& xvelocity,
	SpatialDistribution const& yvelocity,
	EquationOfState const& eos,
	RiemannSolver const& rs,
	PointMotion& pointmotion,
	SourceTerm& external_force,
	OuterBoundary const& obc,
	HydroBoundaryConditions const& hbc,
	bool EntropyCalc,bool CMvalue):
_tessellation(tessellation),
#ifdef RICH_MPI
	_proctess(proctess),
#endif
	_cells(vector<Primitive>()),
	_conservedextensive(vector<Conserved>()),
	_eos(eos),
	_rs(rs),
	_interpolation(interpolation),
	_pointmotion(pointmotion),
	_hbc(hbc),_obc(obc),
	external_force_(external_force),
	_cfl(1./3.),
	_time(0),
	_endtime(-1),
	cycle_(0),
	tracer_(vector<vector<double> >()),
	tracer_flag_(false),
	coldflows_flag_(false),
	densityfloor_(false),
        cfp_(0,0),
	densityMin_(0),
	pressureMin_(0),
	EntropyReCalc_(EntropyCalc),
	_dt_external(0),
  default_pg_(),
  pg_(&default_pg_),
#ifdef RICH_MPI
	procupdate_(0),
#endif
	custom_evolution_manager(),
	custom_evolution_indices(points.size(),0)
{
#ifdef RICH_MPI
	assert(get_mpi_size()%2==0 && "RICH only works with an even number of processes");
#endif

#ifdef RICH_MPI
	_tessellation.Initialise(points,_proctess,&obc);
#else
	_tessellation.Initialise(points, &obc);
#endif
	_cells = InitialiseCells(density, pressure,xvelocity, yvelocity,
		eos, tessellation,CMvalue);

	_conservedextensive = CalcConservedExtensive
	  (CalcConservedIntensive(_cells),tessellation,*pg_);
	_dt_external=-1;
}

hdsim::hdsim(ResetDump const& dump,Tessellation& tessellation,
#ifdef RICH_MPI
	Tessellation& tproc,
#endif
	SpatialReconstruction& interpolation,
	EquationOfState const& eos,RiemannSolver const& rs,
	PointMotion& pointmotion,SourceTerm& external_force,
	OuterBoundary const& obc,HydroBoundaryConditions const& hbc,
	bool EntropyCalc):
_tessellation(tessellation),
#ifdef RICH_MPI
	_proctess(tproc),
#endif
	_cells(dump.snapshot.cells),
	_conservedextensive(vector<Conserved>()),
	_eos(eos),
	_rs(rs),
	_interpolation(interpolation),
	_pointmotion(pointmotion),
	_hbc(hbc),
	_obc(obc),
	external_force_(external_force),
	_cfl(dump.cfl),
	_time(dump.time),
	_endtime(-1),
	cycle_(dump.cycle),
	tracer_(dump.tracers),
	tracer_flag_(!tracer_.empty()),
	coldflows_flag_(dump.coldflows),
	densityfloor_(dump.densityfloor),
        cfp_(dump.a,dump.b),
	densityMin_(dump.densitymin),
	pressureMin_(dump.pressuremin),
	EntropyReCalc_(EntropyCalc),
	_dt_external(-1),
  default_pg_(),
  pg_(&default_pg_),
#ifdef RICH_MPI
	procupdate_(0),
#endif
	custom_evolution_manager(),
	custom_evolution_indices(dump.cevolve)
{
#ifdef RICH_MPI
	assert(get_mpi_size()%2==0 &&
		"RICH only works with an even number of processes");
#endif

#ifndef RICH_MPI
	_tessellation.Initialise(dump.snapshot.mesh_points,&_obc);
#else
	_tessellation.Initialise(dump.snapshot.mesh_points,_proctess,&_obc);
#endif
	_conservedextensive = CalcConservedExtensive
	  (CalcConservedIntensive(_cells),tessellation,*pg_);
}

hdsim::~hdsim(void) {}

void hdsim::SetCfl(double cfl_new)
{
	_cfl = cfl_new;
}

void hdsim::SetTimeStepExternal(double dt)
{
	_dt_external=dt;
}

namespace
{
  vector<char> calc_shocked_cells(const Tessellation& tess,
				  const vector<Primitive>& cells,
				  const HydroBoundaryConditions& hbc,
				  double time)
  {
    vector<char> res(static_cast<size_t>(tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i)
      res[i] = IsShockedCell(tess,static_cast<int>(i),cells,hbc,time) ? 1 : 0;
    return res;
  }

    #ifdef RICH_MPI
    void herd_shocked_status(const Tessellation& tess,
			     vector<char>& shocked_cells)
    {
      vector<char> buf;
      SendRecvShockedStatus(shocked_cells, tess.GetSentPoints(),
			    tess.GetSentProcs(),buf);
      insert_all_to_back(shocked_cells,buf);
    }

  void herd_vector_double(const Tessellation& tess,
			  vector<double>& values)
  {
    vector<double> buf;
    SendRecvVectorDouble(values,tess.GetSentPoints(),tess.GetSentProcs(),buf);
    values = VectorValues(values,tess.GetSelfPoint());
    insert_all_to_back(values,buf);
  }
    #endif

  class ColdFlows
  {
  public:

    ColdFlows(bool active,
	      const Tessellation& tess,
	      const vector<Primitive>& cells,
	      const HydroBoundaryConditions& hbc,
	      double time,
	      const vector<CustomEvolution*>& ce,
	      const vector<double>& g):
      active_(active),
      shocked_cells_(active ? 
		     calc_shocked_cells(tess,cells,hbc,time) : 
		     vector<char>()),
      kinetic_energy_(active ?
		      GetMaxKineticEnergy(tess,cells,ce) :
		      vector<double>()),
      potential_energy_(active ?
			GetForceEnergy(tess,g) : 
			vector<double>()) {}

    #ifdef RICH_MPI
    void herd_data(const Tessellation& tess)
    {
      if(active_){
	herd_shocked_status(tess,shocked_cells_);
	herd_vector_double(tess,kinetic_energy_);
	herd_vector_double(tess,potential_energy_);
      }
    }
    #endif

    /*
    const vector<char>& getShockedStatus(void)
    {
      return shocked_cells_;
    }
    */

    /*
    const vector<double>& getKineticEnergy(void)
    {
      return kinetic_energy_;
    }
    */

    /*
    const vector<double>& getPotentialEnergy(void)
    {
      return potential_energy_;
      }*/

    void fixPressure(vector<Conserved>& intensive,
		     vector<vector<double> >& tracers,
		     const EquationOfState& eos,
		     double as,
		     double bs,
		     vector<CustomEvolution*>& ce,
		     Tessellation& tess,
		     vector<Conserved>& /*extensive*/,
		     bool density_floor)
    {
      if(active_)
	FixPressure(intensive,
		    tracers,
		    eos,
		    kinetic_energy_,
		    potential_energy_,
		    as,
		    bs,
		    ce,
		    tess,
		   // extensive,
		    shocked_cells_,
		    density_floor);
    }

  private:
    const bool active_;
    vector<char> shocked_cells_;
    vector<double> kinetic_energy_;
    vector<double> potential_energy_;
  };

  void substitute_entropy(int num_cells,
			  const vector<Primitive>& cells,
			  const EquationOfState& eos,
			  size_t index,
			  vector<vector<double> >& tracers)
  {
    if(tracers.empty())
      tracers = vector<vector<double> >(static_cast<size_t>(num_cells),vector<double>(1,0));
    for(size_t i=0;i<tracers.size();++i){
      tracers.at(i).at(index) = eos.dp2s(cells.at(i).Density,cells.at(i).Pressure);
    }
  }

  class EdgeLengthCalculator: public Index2Member<double>
  {
  public:

    EdgeLengthCalculator(const Tessellation& tess,
			 const PhysicalGeometry& pg):
      tess_(tess), pg_(pg) {}

    size_t getLength(void) const
    {
      return static_cast<size_t>(tess_.GetTotalSidesNumber());
    }

    double operator()(size_t i) const
    {
      return pg_.calcArea(tess_.GetEdge(static_cast<int>(i)));
    }

  private:
    const Tessellation& tess_;
    const PhysicalGeometry& pg_;
  };
}

void hdsim::TimeAdvance(void)
{
#ifndef RICH_MPI
	PeriodicUpdateCells(_cells,tracer_,custom_evolution_indices,
		_tessellation.GetDuplicatedPoints(),_tessellation.GetTotalPointNumber());
#else
	SendRecvHydro(_cells,custom_evolution_indices,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_eos,_tessellation.GetGhostIndeces()
		,_tessellation.GetTotalPointNumber());

	SendRecvVelocity(_tessellation.GetAllCM(),_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces()
		,_tessellation.GetTotalPointNumber());
#endif

	vector<CustomEvolution*> custom_evolutions =
		convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	vector<Vector2D> point_velocity = 
	  _pointmotion.calcAllVelocities
	  (_tessellation,_cells,_time,custom_evolutions);

#ifndef RICH_MPI
	PeriodicVelocityExchange(point_velocity,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetTotalPointNumber());
#else
	SendRecvVelocity(point_velocity,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces(),
		_tessellation.GetTotalPointNumber());
#endif

	const vector<Vector2D> fv=
	  _tessellation.calc_edge_velocities
	  (&_hbc,point_velocity,_time);

	const double dt = determine_time_step(_cfl*CalcTimeStep(_tessellation, _cells, fv,_hbc,
								_time,custom_evolutions),
					      _dt_external, _time, _endtime);

	if(coldflows_flag_)
	  substitute_entropy(_tessellation.GetPointNo(),
			     _cells,_eos,0,tracer_);

#ifdef RICH_MPI
	if(tracer_flag_)
		SendRecvTracers(tracer_,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces()
		,_tessellation.GetTotalPointNumber());
#endif
	vector<Conserved> fluxes = calc_fluxes
	  (_tessellation,  _cells, dt, _time,
	   _interpolation,
	   fv, _hbc, _rs,
	   custom_evolutions,
	   custom_evolution_manager,
	   tracer_);					

	const vector<double> lengths = serial_generate(EdgeLengthCalculator(_tessellation,*pg_));

	vector<vector<double> > tracer_extensive = 
	  calc_extensive_tracer(tracer_,
				_tessellation,
				_cells,
				*pg_);
 
	really_update_extensive_tracers(tracer_extensive,
					tracer_,
					_cells,
					_tessellation,
					fluxes,
					_time,dt,
					_hbc,
					_interpolation,
					custom_evolutions,
					custom_evolution_manager,
					fv,lengths);

	vector<double> g;
	ExternalForceContribution(_tessellation,*pg_,
				  _cells,external_force_,_time,dt,
		_conservedextensive,_hbc,fluxes,point_velocity,g,coldflows_flag_,tracer_,
		lengths);

	ColdFlows cold_flows(coldflows_flag_,
			     _tessellation,
			     _cells,
			     _hbc,
			     _time,
			     custom_evolutions,
			     g);

	UpdateConservedExtensive(_tessellation, fluxes, dt,
		_conservedextensive,_hbc,lengths);

#ifndef RICH_MPI
	MoveMeshPoints(point_velocity, dt, _tessellation);
#else
	if(procupdate_)
		procupdate_->Update(_proctess,_tessellation);
	MoveMeshPoints(point_velocity, dt, _tessellation,_proctess);
#endif

#ifdef RICH_MPI
	vector<Conserved> ptoadd;
	vector<vector<double> > ttoadd;
	vector<size_t> ctoadd;
	SendRecvExtensive(_conservedextensive,tracer_extensive,custom_evolution_indices,
		_tessellation.GetSentPoints(),_tessellation.GetSentProcs(),ptoadd,ttoadd,
		ctoadd);

	KeepLocalPoints(_conservedextensive,tracer_extensive,custom_evolution_indices,
		_tessellation.GetSelfPoint());

	insert_all_to_back(_conservedextensive,ptoadd);
	insert_all_to_back(tracer_extensive,ttoadd);
	insert_all_to_back(custom_evolution_indices,ctoadd);

	custom_evolutions =
		convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	cold_flows.herd_data(_tessellation);
#endif

	vector<Conserved> intensive = 
	  calc_conserved_intensive(_tessellation,
				   _conservedextensive,
				   *pg_);

	cold_flows.fixPressure(intensive,
			       tracer_extensive,
			       _eos,
			       cfp_.as, cfp_.bs,
			       custom_evolutions,
			       _tessellation,
			       _conservedextensive,
			       densityfloor_);

	UpdatePrimitives(intensive, 
			 _eos, _cells,custom_evolutions,_cells,
		densityfloor_,densityMin_,pressureMin_,_tessellation,_time,
		tracer_extensive);

	MakeTracerIntensive(tracer_,tracer_extensive,
			    _tessellation,_cells, *pg_);

	_time += dt;
	cycle_++;
}

namespace {

  template<class T> class AverageCalculator: public Index2Member<T>
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

void hdsim::TimeAdvanceElad(void)
{
#ifndef RICH_MPI
	PeriodicUpdateCells(_cells,tracer_,custom_evolution_indices,
		_tessellation.GetDuplicatedPoints(),_tessellation.GetTotalPointNumber());
#else
	SendRecvHydro(_cells,custom_evolution_indices,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_eos,_tessellation.GetGhostIndeces()
		,_tessellation.GetTotalPointNumber());

	SendRecvVelocity(_tessellation.GetAllCM(),_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces()
		,_tessellation.GetTotalPointNumber());
#endif

	vector<CustomEvolution*> custom_evolutions =
		convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	vector<Vector2D> point_velocity = 
	  _pointmotion.calcAllVelocities
	  (_tessellation,_cells,_time,custom_evolutions);

#ifndef RICH_MPI
	PeriodicVelocityExchange(point_velocity,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetTotalPointNumber());
#else
	SendRecvVelocity(point_velocity,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces(),
		_tessellation.GetTotalPointNumber());
#endif

	const vector<Vector2D> fv0=
	  _tessellation.calc_edge_velocities
	  (&_hbc,point_velocity,_time);

	const double dt = determine_time_step(_cfl*CalcTimeStep(_tessellation, _cells, fv0,_hbc,
								_time,custom_evolutions),
					      _dt_external, _time, _endtime);

	if(coldflows_flag_)
	  substitute_entropy(_tessellation.GetPointNo(),
			     _cells,_eos,0,tracer_);

#ifdef RICH_MPI
	if(tracer_flag_)
		SendRecvTracers(tracer_,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces()
		,_tessellation.GetTotalPointNumber());
#endif
	vector<Conserved> fluxes_0 = calc_fluxes
	  (_tessellation,  _cells, dt, _time,
	   _interpolation,
	   fv0, _hbc, _rs,
	   custom_evolutions,
	   custom_evolution_manager,
	   tracer_);

	vector<Conserved> fluxes_1 = calc_fluxes
	  (_tessellation,  _cells, dt, _time,
	   _interpolation,
	   serial_generate
	   (FaceVertexVelocityCalculator
	    (_tessellation,
	     point_velocity,
	     &std::pair<Vector2D,Vector2D>::first,
	     fv0,_hbc)),
	   _hbc, _rs,
	   custom_evolutions,
	   custom_evolution_manager,
	   tracer_);

	vector<Conserved> fluxes_2 = calc_fluxes
	  (_tessellation,  _cells, dt, _time,
	   _interpolation,
	   serial_generate
	   (FaceVertexVelocityCalculator
	    (_tessellation,
	     point_velocity,
	     &std::pair<Vector2D,Vector2D>::second,
	     fv0,_hbc)),
	   _hbc, _rs,
	   custom_evolutions,
	   custom_evolution_manager,
	   tracer_);

	vector<Conserved> fluxes = 
	  average(fluxes_0,average(fluxes_1,fluxes_2));

	const vector<double> lengths = serial_generate(EdgeLengthCalculator(_tessellation,*pg_));

	vector<vector<double> > tracer_extensive = 
	  calc_extensive_tracer(tracer_,
				_tessellation,
				_cells,
				*pg_);
 
	really_update_extensive_tracers(tracer_extensive,
					tracer_,
					_cells,
					_tessellation,
					fluxes,
					_time,dt,
					_hbc,
					_interpolation,
					custom_evolutions,
					custom_evolution_manager,
					fv0,lengths);

	vector<double> g;
	ExternalForceContribution(_tessellation,*pg_,
				  _cells,external_force_,_time,dt,
		_conservedextensive,_hbc,fluxes,point_velocity,g,coldflows_flag_,tracer_,
		lengths);

	ColdFlows cold_flows(coldflows_flag_,
			     _tessellation,
			     _cells,
			     _hbc,
			     _time,
			     custom_evolutions,
			     g);

	UpdateConservedExtensive(_tessellation, fluxes, dt,
		_conservedextensive,_hbc,lengths);

#ifndef RICH_MPI
	MoveMeshPoints(point_velocity, dt, _tessellation);
#else
	if(procupdate_)
		procupdate_->Update(_proctess,_tessellation);
	MoveMeshPoints(point_velocity, dt, _tessellation,_proctess);
#endif

#ifdef RICH_MPI
	vector<Conserved> ptoadd;
	vector<vector<double> > ttoadd;
	vector<size_t> ctoadd;
	SendRecvExtensive(_conservedextensive,tracer_extensive,custom_evolution_indices,
		_tessellation.GetSentPoints(),_tessellation.GetSentProcs(),ptoadd,ttoadd,
		ctoadd);

	KeepLocalPoints(_conservedextensive,tracer_extensive,custom_evolution_indices,
		_tessellation.GetSelfPoint());

	insert_all_to_back(_conservedextensive,ptoadd);
	insert_all_to_back(tracer_extensive,ttoadd);
	insert_all_to_back(custom_evolution_indices,ctoadd);

	custom_evolutions =
		convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	cold_flows.herd_data(_tessellation);
#endif

	vector<Conserved> intensive = 
	  calc_conserved_intensive(_tessellation,
				   _conservedextensive,
				   *pg_);

	cold_flows.fixPressure(intensive,
			       tracer_extensive,
			       _eos,
			       cfp_.as, cfp_.bs,
			       custom_evolutions,
			       _tessellation,
			       _conservedextensive,
			       densityfloor_);

	UpdatePrimitives(intensive, 
			 _eos, _cells,custom_evolutions,_cells,
		densityfloor_,densityMin_,pressureMin_,_tessellation,_time,
		tracer_extensive);

	MakeTracerIntensive(tracer_,tracer_extensive,
			    _tessellation,_cells, *pg_);

	_time += dt;
	cycle_++;
}

const PhysicalGeometry& hdsim::getPhysicalGeometry(void) const
{
  return *pg_;
}

const Tessellation& hdsim::GetTessellation(void) const
{
	return _tessellation;
}

#ifdef RICH_MPI
Tessellation const& hdsim::GetProcTessellation(void) const
{
	return _proctess;
}
#endif

void hdsim::changePhysicalGeometry(const PhysicalGeometry* pg)
{
  pg_ = pg;
  _conservedextensive = CalcConservedExtensive
    (CalcConservedIntensive(_cells),_tessellation,*pg_);
}

void hdsim::TimeAdvance2Mid(void)
{
	const double dt=TimeAdvance2mid
		(_tessellation,
#ifdef RICH_MPI
		_proctess,
#endif
		_cells,_pointmotion,
		_hbc,_interpolation,_rs,_eos,external_force_,_time,_cfl,_endtime,
		tracer_,_dt_external,custom_evolution_indices,
		 custom_evolution_manager, *pg_,
#ifdef RICH_MPI
		procupdate_,
#endif
		tracer_flag_,
		coldflows_flag_,cfp_.as,cfp_.bs,densityfloor_,densityMin_,pressureMin_,
		EntropyReCalc_);

	_time += dt;
	++cycle_;
}

void hdsim::TimeAdvanceElad2(void)
{

	// create ghost points if needed
#ifndef RICH_MPI
	PeriodicUpdateCells(_cells,tracer_,custom_evolution_indices,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetTotalPointNumber());
#else
	SendRecvHydro(_cells,custom_evolution_indices,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_eos,_tessellation.GetGhostIndeces(),_tessellation.GetTotalPointNumber());

	SendRecvVelocity(_tessellation.GetAllCM(),_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces(),_tessellation.GetTotalPointNumber());
#endif
	vector<CustomEvolution*> CellsEvolve=convert_indices_to_custom_evolution(
		custom_evolution_manager,custom_evolution_indices);
	vector<Vector2D> oldpoints=_tessellation.GetMeshPoints();
	oldpoints.resize(static_cast<size_t>(_tessellation.GetPointNo()));

	//do half time step
	vector<Vector2D> point_velocities = calc_point_velocities
		(_tessellation, _cells, _pointmotion, _time,CellsEvolve);

#ifndef RICH_MPI
	PeriodicVelocityExchange(point_velocities,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetTotalPointNumber());
#else
	SendRecvVelocity(point_velocities,_tessellation.GetDuplicatedPoints(),
			 _tessellation.GetDuplicatedProcs(),
			 _tessellation.GetGhostIndeces(),
			 _tessellation.GetTotalPointNumber());
#endif

	vector<Vector2D> edge_velocities =_tessellation.calc_edge_velocities
		(&_hbc,point_velocities,_time);

	double dt = determine_time_step
		(_cfl*CalcTimeStep(_tessellation,_cells,edge_velocities,_hbc,_time,CellsEvolve),
		 _dt_external,_time,_endtime);

	// Entropy and tracers evolution
	vector<double> g,Ek,Ef;
	vector<char> shockedcells;
	if(coldflows_flag_)
	{
		const int n=_tessellation.GetPointNo();
		for(int i=0;i<n;++i)
			tracer_[static_cast<size_t>(i)][0]=_eos.dp2s(_cells[static_cast<size_t>(i)].Density,
							_cells[static_cast<size_t>(i)].Pressure);
#ifdef RICH_MPI
		SendRecvTracers(tracer_,_tessellation.GetDuplicatedPoints(),
			_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces(),_tessellation.GetTotalPointNumber());
		if(tracer_.size()!=_cells.size())
		{
			UniversalError eo("Tracers and cells have different length at first half time step");
			eo.AddEntry("CPU rank",get_mpi_rank());
			throw eo;
		}
#endif
		shockedcells.resize(static_cast<size_t>(n));
		for(int i=0;i<n;++i)
			shockedcells[static_cast<size_t>(i)]=IsShockedCell(_tessellation,i,_cells,_hbc,_time) ? 1 : 0 ;
	}
#ifdef RICH_MPI
	if(tracer_flag_&&!coldflows_flag_)
	{
		SendRecvTracers(tracer_,_tessellation.GetDuplicatedPoints(),
			_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces(),_tessellation.GetTotalPointNumber());
		if(tracer_.size()!=_cells.size())
		{
			UniversalError eo("Tracers and cells have different length at first half time step");
			eo.AddEntry("CPU rank",get_mpi_rank());
			throw eo;
		}
	}
#endif

	vector<Conserved> fluxes0 = calc_fluxes
		(_tessellation, _cells, 0.5*dt, _time, _interpolation,
		edge_velocities, _hbc, _rs,CellsEvolve,
		custom_evolution_manager,tracer_);

	vector<Conserved> fluxes1 = calc_fluxes
	  (_tessellation, _cells, 0.5*dt, _time, _interpolation,
	   serial_generate
	   (FaceVertexVelocityCalculator
	    (_tessellation,
	     point_velocities,
	     &std::pair<Vector2D,Vector2D>::first,
	     edge_velocities,_hbc)),
	   _hbc, _rs,CellsEvolve,
	   custom_evolution_manager,tracer_);

	vector<Conserved> fluxes2 = calc_fluxes
	  (_tessellation, _cells, 0.5*dt, _time, _interpolation,
	   serial_generate
	   (FaceVertexVelocityCalculator
	    (_tessellation,
	     point_velocities,
	     &std::pair<Vector2D,Vector2D>::second,
	     edge_velocities,_hbc)),
	   _hbc, _rs,CellsEvolve,
	   custom_evolution_manager,tracer_);

	vector<Conserved> fluxes =
	  average(fluxes0,average(fluxes1,fluxes2));
	
	vector<double> lengths;
	int nsides=_tessellation.GetTotalSidesNumber();
	lengths.resize(static_cast<size_t>(nsides));
	for(int i=0;i<nsides;++i)
		lengths[static_cast<size_t>(i)]=_tessellation.GetEdge(i).GetLength();

	vector<vector<double> > old_trace;
	vector<vector<double> > tracer_extensive;
	if(tracer_flag_)
	{
		vector<vector<double> > trace_change;
		trace_change = CalcTraceChange
			(tracer_,_cells,_tessellation,fluxes,0.5*dt,_hbc,
			 _interpolation,_time,CellsEvolve,
			custom_evolution_manager,
			edge_velocities,lengths);
		MakeTracerExtensive(tracer_,_tessellation,_cells,tracer_extensive);
		old_trace=tracer_extensive;
		UpdateTracerExtensive(tracer_extensive,trace_change,CellsEvolve,_cells,
			_tessellation,_time);
	}

	vector<Conserved> intensive = CalcConservedIntensive(_cells);

	vector<Conserved> extensive = CalcConservedExtensive
		(intensive, _tessellation, *pg_);

	// Save extensive variables of beginning of time step
	// if(traceflag)
	// MakeTracerExtensive(tracers,tess,cells,old_trace);
	vector<Conserved> old_extensive=extensive;

	UpdateConservedExtensive(_tessellation, fluxes, 0.5*dt,
		extensive, _hbc,lengths);

	ExternalForceContribution
		(_tessellation, *pg_, _cells,external_force_, _time, 0.5*dt,
		extensive, _hbc,fluxes,point_velocities,g,coldflows_flag_,tracer_,lengths);

	if(coldflows_flag_)
	{
		Ek=GetMaxKineticEnergy(_tessellation,_cells,CellsEvolve);
		Ef=GetForceEnergy(_tessellation,g);
	}

#ifndef RICH_MPI
	MoveMeshPoints(point_velocities, 0.5*dt, _tessellation);
#else
	//	if(procupdate!=0)
	//		procupdate->Update(proctess,tess);
	MoveMeshPoints(point_velocities,0.5*dt, _tessellation,_proctess);
#endif

#ifdef RICH_MPI
	vector<Conserved> ptoadd;
	vector<vector<double> > ttoadd;
	vector<size_t> ctemp(custom_evolution_indices);
	vector<size_t> ctoadd;
	SendRecvExtensive(extensive,tracer_extensive,custom_evolution_indices,
		_tessellation.GetSentPoints(),_tessellation.GetSentProcs(),ptoadd,ttoadd,
		ctoadd);

	KeepLocalPoints(extensive,tracer_extensive,custom_evolution_indices,
		_tessellation.GetSelfPoint());

	if(!ptoadd.empty())
	{
		extensive.insert(extensive.end(),ptoadd.begin(),
			ptoadd.end());
	}

	if(!ttoadd.empty())
	{
		tracer_extensive.insert(tracer_extensive.end(),ttoadd.begin(),ttoadd.end());
	}

	if(!ctoadd.empty())
	{
		custom_evolution_indices.insert(custom_evolution_indices.end(),ctoadd.begin(),
			ctoadd.end());
	}

	SendRecvExtensive(old_extensive,old_trace,ctemp,
		_tessellation.GetSentPoints(),_tessellation.GetSentProcs(),ptoadd,ttoadd,
		ctoadd);

	KeepLocalPoints(old_extensive,old_trace,ctemp,_tessellation.GetSelfPoint());

	if(!ptoadd.empty())
	{
		old_extensive.insert(old_extensive.end(),ptoadd.begin(),
			ptoadd.end());
	}

	if(!ttoadd.empty())
	{
		old_trace.insert(old_trace.end(),ttoadd.begin(),ttoadd.end());
	}

	vector<Vector2D> vtoadd;
	SendRecvOldVector2D(oldpoints,_tessellation.GetSentPoints(),
		_tessellation.GetSentProcs(),vtoadd);
	oldpoints=VectorValues(oldpoints,_tessellation.GetSelfPoint());
	if(!vtoadd.empty())
		oldpoints.insert(oldpoints.end(),vtoadd.begin(),vtoadd.end());

	CellsEvolve =
		convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	if(coldflows_flag_)
	{
		vector<char> btoadd;
		SendRecvShockedStatus(shockedcells,_tessellation.GetSentPoints(),_tessellation.GetSentProcs(),
			btoadd);
		shockedcells=VectorValues(shockedcells,_tessellation.GetSelfPoint());
		if(!btoadd.empty())
			shockedcells.insert(shockedcells.end(),btoadd.begin(),btoadd.end());
		vector<double> Ekadd;
		SendRecvVectorDouble(Ek,_tessellation.GetSentPoints(),_tessellation.GetSentProcs(),
			Ekadd);
		Ek=VectorValues(Ek,_tessellation.GetSelfPoint());
		if(!Ekadd.empty())
			Ek.insert(Ek.end(),Ekadd.begin(),Ekadd.end());
		SendRecvVectorDouble(Ef,_tessellation.GetSentPoints(),_tessellation.GetSentProcs(),
			Ekadd);
		Ef=VectorValues(Ef,_tessellation.GetSelfPoint());
		if(!Ekadd.empty())
			Ef.insert(Ef.end(),Ekadd.begin(),Ekadd.end());
	}
#endif

	UpdateConservedIntensive(_tessellation, extensive, intensive);

	if(coldflows_flag_)
		FixPressure(intensive,tracer_extensive,_eos,Ek,Ef,
			    cfp_.as,cfp_.bs,CellsEvolve,
			    _tessellation,/*extensive,*/shockedcells,densityfloor_);

	_cells.resize(static_cast<size_t>(_tessellation.GetPointNo()));
	UpdatePrimitives(intensive, _eos, _cells,CellsEvolve,_cells,densityfloor_,
			 densityMin_,pressureMin_,_tessellation,_time+0.5*dt,tracer_);
	if(tracer_flag_)
	{
		MakeTracerIntensive(tracer_,tracer_extensive,_tessellation,_cells,*pg_);
	}

	// End half step

	// create ghost points if needed
#ifndef RICH_MPI
	PeriodicUpdateCells(_cells,tracer_,custom_evolution_indices,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetTotalPointNumber());
#else
	SendRecvHydro(_cells,custom_evolution_indices,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_eos,_tessellation.GetGhostIndeces(),_tessellation.GetTotalPointNumber());

	SendRecvVelocity(_tessellation.GetAllCM(),_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces(),_tessellation.GetTotalPointNumber());
#endif

	CellsEvolve=convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	point_velocities = calc_point_velocities(_tessellation,_cells,
		_pointmotion, _time+0.5*dt,CellsEvolve);

#ifndef RICH_MPI
	PeriodicVelocityExchange(point_velocities,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetTotalPointNumber());
#else
	SendRecvVelocity(point_velocities,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces(),_tessellation.GetTotalPointNumber());
#endif

	edge_velocities = _tessellation.calc_edge_velocities(&_hbc,point_velocities,_time);

	if(coldflows_flag_&&EntropyReCalc_)
	{
		const int n=_tessellation.GetPointNo();
		for(size_t i=0;i<static_cast<size_t>(n);++i)
			tracer_[i][0]=_eos.dp2s(_cells[i].Density,_cells[i].Pressure);
#ifdef RICH_MPI
		SendRecvTracers(tracer_,_tessellation.GetDuplicatedPoints(),
			_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces(),_tessellation.GetTotalPointNumber());
		if(tracer_.size()!=_cells.size())
		{
			UniversalError eo("Tracers and cells have different length at second half time step");
			eo.AddEntry("CPU rank",get_mpi_rank());
			throw eo;
		}

#endif
	}

#ifdef RICH_MPI
	if(tracer_flag_&&(!coldflows_flag_||!EntropyReCalc_))
	{
		SendRecvTracers(tracer_,_tessellation.GetDuplicatedPoints(),
			_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces(),_tessellation.GetTotalPointNumber());
		if(tracer_.size()!=_cells.size())
		{
			UniversalError eo("Tracers and cells have different length at second half time step");
			eo.AddEntry("CPU rank",get_mpi_rank());
			throw eo;
		}
	}
#endif
	fluxes0 = calc_fluxes
		(_tessellation, _cells, dt, _time+dt/2, _interpolation,
		edge_velocities, _hbc, _rs,CellsEvolve,
		custom_evolution_manager,tracer_);

	fluxes1 = calc_fluxes
	  (_tessellation, _cells, dt, _time+dt/2, _interpolation,
	   serial_generate
	   (FaceVertexVelocityCalculator
	    (_tessellation,
	     point_velocities,
	     &std::pair<Vector2D,Vector2D>::first,
	     edge_velocities,_hbc)),
	   _hbc, _rs,CellsEvolve,
	   custom_evolution_manager,tracer_);

	fluxes2 = calc_fluxes
	  (_tessellation, _cells, dt, _time+dt/2, _interpolation,
	   serial_generate
	   (FaceVertexVelocityCalculator
	    (_tessellation,
	     point_velocities,
	     &std::pair<Vector2D,Vector2D>::second,
	     edge_velocities,_hbc)),
	   _hbc, _rs,CellsEvolve,
	   custom_evolution_manager,tracer_);

	fluxes = average(fluxes0,
			 average(fluxes1,
				 fluxes2));
	
	if(coldflows_flag_)
	{
		const int n=_tessellation.GetPointNo();
		shockedcells.resize(static_cast<size_t>(n));
		for(int i=0;i<n;++i)
			shockedcells[static_cast<size_t>(i)]=IsShockedCell(_tessellation,i,_cells,_hbc,_time) ? 1 : 0 ;
	}

	nsides=_tessellation.GetTotalSidesNumber();
	lengths.resize(static_cast<size_t>(nsides));
	for(int i=0;i<nsides;++i)
		lengths[static_cast<size_t>(i)]=_tessellation.GetEdge(i).GetLength();

	if(tracer_flag_)
	{
		vector<vector<double> > trace_change;
		trace_change = CalcTraceChange
		  (tracer_,_cells,_tessellation,fluxes,dt,_hbc,
			_interpolation,_time,CellsEvolve,
			custom_evolution_manager,
			edge_velocities,lengths);
		UpdateTracerExtensive(old_trace,trace_change,CellsEvolve,_cells,
			_tessellation,_time);
	}

	UpdateConservedExtensive(_tessellation, fluxes, dt,
		old_extensive, _hbc,lengths);

	ExternalForceContribution
		(_tessellation, *pg_, _cells,external_force_, _time+0.5*dt,dt,
		old_extensive, _hbc,fluxes,point_velocities,g,coldflows_flag_,tracer_,lengths);

	if(coldflows_flag_)
	{
		Ek=GetMaxKineticEnergy(_tessellation,_cells,CellsEvolve);
		Ef=GetForceEnergy(_tessellation,g);
	}
#ifndef RICH_MPI
	MoveMeshPoints(point_velocities, dt, _tessellation,oldpoints);
#else
	if(procupdate_!=0)
		procupdate_->Update(_proctess,_tessellation);
	MoveMeshPoints(point_velocities, dt, _tessellation,_proctess,oldpoints);
#endif

#ifdef RICH_MPI
	SendRecvExtensive(old_extensive,old_trace,custom_evolution_indices,
		_tessellation.GetSentPoints(),_tessellation.GetSentProcs(),ptoadd,ttoadd,
		ctoadd);

	KeepLocalPoints(old_extensive,old_trace,custom_evolution_indices,
		_tessellation.GetSelfPoint());

	if(!ptoadd.empty())
	{
		old_extensive.insert(old_extensive.end(),ptoadd.begin(),
			ptoadd.end());
	}

	if(!ttoadd.empty())
	{
		old_trace.insert(old_trace.end(),ttoadd.begin(),ttoadd.end());
	}

	if(!ctoadd.empty())
	{
		custom_evolution_indices.insert(custom_evolution_indices.end(),ctoadd.begin(),
			ctoadd.end());
	}

	CellsEvolve =
		convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	if(coldflows_flag_)
	{
		vector<char> btoadd;
		SendRecvShockedStatus(shockedcells,_tessellation.GetSentPoints(),_tessellation.GetSentProcs(),
			btoadd);
		shockedcells=VectorValues(shockedcells,_tessellation.GetSelfPoint());
		if(!btoadd.empty())
			shockedcells.insert(shockedcells.end(),btoadd.begin(),btoadd.end());
		vector<double> Ekadd;
		SendRecvVectorDouble(Ek,_tessellation.GetSentPoints(),_tessellation.GetSentProcs(),
			Ekadd);
		Ek=VectorValues(Ek,_tessellation.GetSelfPoint());
		if(!Ekadd.empty())
			Ek.insert(Ek.end(),Ekadd.begin(),Ekadd.end());
		SendRecvVectorDouble(Ef,_tessellation.GetSentPoints(),_tessellation.GetSentProcs(),
			Ekadd);
		Ef=VectorValues(Ef,_tessellation.GetSelfPoint());
		if(!Ekadd.empty())
			Ef.insert(Ef.end(),Ekadd.begin(),Ekadd.end());
	}
#endif

	UpdateConservedIntensive(_tessellation, old_extensive, intensive);

	if(coldflows_flag_)
		FixPressure(intensive,old_trace,_eos,Ek,Ef,cfp_.as,cfp_.bs,CellsEvolve,_tessellation,
		/*extensive,*/shockedcells,densityfloor_);

	UpdatePrimitives
		(intensive, _eos, _cells,CellsEvolve,_cells,densityfloor_,
		 densityMin_,pressureMin_,_tessellation,_time+dt,tracer_);

	if(tracer_flag_)
	{
		MakeTracerIntensive(tracer_,old_trace,_tessellation,_cells,*pg_);
	}
	
	_time += dt;
	++cycle_;
}

void hdsim::addTracer(SpatialDistribution const& tp)
{
	const int n = _tessellation.GetPointNo();
	if(!tracer_flag_){
	  tracer_.resize(static_cast<size_t>(n));
		tracer_flag_ = true;
	}

	for(int i=0;i<n;++i)
	  tracer_[static_cast<size_t>(i)].push_back(tp(_tessellation.GetCellCM(i)));
}

// Diagnostics

int hdsim::GetEdgeNo(void) const
{
	return _tessellation.GetTotalSidesNumber();
}

Edge hdsim::GetEdge(int i) const
{
	return _tessellation.GetEdge(i);
}

int hdsim::GetCellNo(void) const
{
	return _tessellation.GetPointNo();
}

Primitive hdsim::GetCell(int i) const
{
  assert(i>=0);
  return _cells.at(static_cast<size_t>(i));
}

Vector2D hdsim::GetMeshPoint(int i) const
{
	return _tessellation.GetMeshPoint(i);
}

double hdsim::GetTime(void) const
{
	return _time;
}

double hdsim::GetCellVolume(int index) const
{
	return _tessellation.GetVolume(index);
}

int hdsim::GetCycle(void) const
{
	return cycle_;
}

void hdsim::SetEndTime(double endtime)
{
	_endtime=endtime;
}

void hdsim::SetColdFlows(double as,double bs)
{
	const int n = _tessellation.GetPointNo();

	coldflows_flag_=true;
	if(!tracer_flag_)
	{
	  tracer_.resize(static_cast<size_t>(n));
		tracer_flag_ = true;
	}

	for(int i=0;i<n;++i)
	  tracer_[static_cast<size_t>(i)].push_back(0);
	cfp_.as=as;
	cfp_.bs=bs;
}

vector<vector<double> > const& hdsim::getTracers(void) const
{
	return tracer_;
}

vector<vector<double> >& hdsim::getTracers(void)
{
	return tracer_;
}

void hdsim::TracerReset(double alpha,SpatialDistribution const& originalD,
	SpatialDistribution const& originalP,SpatialDistribution const& originalVx,
	SpatialDistribution const& originalVy,vector<SpatialDistribution const*> const& originalTracers, int tracerindex)
{
	vector<CustomEvolution*> cevolve=convert_indices_to_custom_evolution(
		custom_evolution_manager,custom_evolution_indices);
	TracerResetCalc(alpha,originalD,originalP,originalVx,originalVy,originalTracers,
		_cells,_tessellation,tracer_,tracerindex,_eos,cevolve,coldflows_flag_);
	return;
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
				if(rank<procorder[i])
				{
				  if(toremove[static_cast<size_t>(index)].empty())
					  MPI_Send(&temp,1,MPI_INT,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD);
					else
					  MPI_Send(&toremove[index][0],static_cast<int>(toremove[static_cast<size_t>(index)].size()),
						   MPI_INT,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD);
					MPI_Probe(procorder[static_cast<size_t>(i)],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					if(status.MPI_TAG==1)
						MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					else
					{
						int count;
						MPI_Get_count(&status,MPI_INT,&count);
						recv[index].resize(count);
						MPI_Recv(&recv[index][0],count,MPI_INT,procorder[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					}
				}
				else
				{
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					if(status.MPI_TAG==1)
						MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					else
					{
						int count;
						MPI_Get_count(&status,MPI_INT,&count);
						recv[index].resize(count);
						MPI_Recv(&recv[index][0],count,MPI_INT,procorder[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					}
					if(toremove[index].empty())
						MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
					else
					  MPI_Send(&toremove[index][0],static_cast<int>(toremove[index].size()),
						MPI_INT,procorder[i],0,MPI_COMM_WORLD);
				}
			}
		}
		for(int i=0;i<nlist;++i)
			if(!recv[i].empty())
				RemoveVector(nghost[i],recv[i]);
	}
#endif
}

namespace {
  class VolumeExtractor: public Index2Member<double>
  {
  public:

    VolumeExtractor(const Tessellation& tess):
      tess_(tess) {}

    size_t getLength(void) const
    {
      return static_cast<size_t>(tess_.GetPointNo());
    }

    double operator()(size_t i) const
    {
      return tess_.GetVolume(static_cast<int>(i));
    }

  private:
    const Tessellation& tess_;
  };

  vector<int> calc_TotalNeigh(const vector<vector<int> >& vol_index)
  {
    vector<int> res;
    for(size_t i=0;i<vol_index.size();++i)
      res.insert(res.end(),
		 vol_index[i].begin(),
		 vol_index[i].end());
    sort(res.begin(),res.end());
    return unique(res);
  }

  vector<Conserved> calc_c_temp(const vector<int>& total_neigh,
				const vector<vector<int> >& vol_index,
				const vector<vector<double> >& dv,
				const vector<Primitive>& cells,
				#ifdef RICH_MPI
				const vector<Primitive>& mpi_cells,
				#endif
				const vector<int>& to_remove,
				int n)
  {
    vector<Conserved> res(total_neigh.size());
    for(size_t i=0;i<dv.size();++i){
      for(size_t j=0;j<vol_index[i].size();++j){
	const size_t index = static_cast<size_t>(lower_bound(total_neigh.begin(),total_neigh.end(),
						  vol_index[i][j])-total_neigh.begin());
	if(i<static_cast<size_t>(n))
	  res[index] += Primitive2Conserved(cells[static_cast<size_t>(to_remove[i])],dv[i][j]);
#ifdef RICH_MPI
	else{
	  res[index] += Primitive2Conserved(mpi_cells[i-static_cast<size_t>(n)],dv[i][j]);
	}
#endif
      }
    }
    return res;
  }

  vector<vector<double> > calc_t_temp(bool traceractive,
				      const vector<int>& total_neigh,
				      const vector<vector<int> >& vol_index,
				      const vector<vector<double> >& dv,
				      const vector<Primitive>& cells,
				      const vector<vector<double> >& tracers,
#ifdef RICH_MPI
				      const vector<Primitive>& mpi_cells,
				      const vector<vector<double> >& mpi_tracer,
#endif
				      const vector<int>& to_remove,
				      int n)
				      
  {
    if(!traceractive)
      return vector<vector<double> >();
    vector<vector<double> > res(total_neigh.size(),vector<double>(tracers[0].size(),0));
    for(size_t i=0;i<dv.size();++i){
      for(size_t j=0;j<vol_index[i].size();++j){
	const size_t index = static_cast<size_t>(lower_bound(total_neigh.begin(),total_neigh.end(),
						  vol_index[i][j])-total_neigh.begin());
				      
	if(i<static_cast<size_t>(n)){
	  for(size_t k=0;k<res[0].size();++k)
	    res[index][k] += dv[i][j]*tracers[static_cast<size_t>(to_remove[i])][k]*cells[static_cast<size_t>(to_remove[i])].Density;
	}
	#ifdef RICH_MPI
	else{
	  for(size_t k=0;k<res[0].size();++k)
	    res[index][k] += dv[i][j]*mpi_tracer[i-static_cast<size_t>(n)][k]*mpi_cells[i-static_cast<size_t>(n)].Density;
	}	  
	#endif
      }
    }
    return res;
  }

  void remove_update_cells_tracers(const bool tracer_active,
				   const vector<int>& total_neigh,
				   const vector<int>& to_remove,
				   const Tessellation& tess,
				   const vector<double>& old_volumes,
				   const vector<Conserved>& c_temp,
				   const vector<vector<double> >& t_temp,
				   const EquationOfState& eos,
				   vector<Primitive>& cells,
				   vector<vector<double> >& tracers)
  {
    for(size_t i=0;i<total_neigh.size();++i){
      const size_t index = static_cast<size_t>(lower_bound(to_remove.begin(),to_remove.end(),
						total_neigh[i])-to_remove.begin());
      const double volume = tess.GetVolume(total_neigh[i]-static_cast<int>(index));
      const Conserved old_extensive = Primitive2Conserved(cells[static_cast<size_t>(total_neigh[i])],
							  old_volumes[static_cast<size_t>(total_neigh[i])]);
      const double old_density = cells[static_cast<size_t>(total_neigh[i])].Density;
      cells[static_cast<size_t>(total_neigh[i])] = Conserved2Primitive((c_temp[i]+old_extensive)/volume,eos);
      if(tracer_active){
	const double new_density = cells[static_cast<size_t>(total_neigh[i])].Density;
	tracers[static_cast<size_t>(total_neigh[i])] = (1./(new_density*volume))*
	  (t_temp[i]+old_volumes[static_cast<size_t>(total_neigh[i])]*
	   old_density*tracers[static_cast<size_t>(total_neigh[i])]);	  
      }
    }
  }
}

vector<int> hdsim::RemoveCells(RemovalStrategy const* remove)
{
	if(!remove)
		throw UniversalError("No Removal strategy");
	const bool traceractive = !tracer_.empty();
	vector<int> ToRemove=remove->CellsToRemove(_tessellation,_cells,tracer_,_time);
	//	int n=int(ToRemove.size());
	if(!ToRemove.empty())
		sort(ToRemove.begin(),ToRemove.end());
	const vector<double> OldVol = serial_generate
	  (VolumeExtractor(_tessellation));
	// Change the tessellation
	vector<vector<int> > VolIndex;
	vector<vector<double> > dv;
	_tessellation.RemoveCells(ToRemove,VolIndex,dv);
	int n=static_cast<int>(lower_bound(ToRemove.begin(),ToRemove.end(),static_cast<int>(OldVol.size()))- ToRemove.begin());
	// gather all the relevant neighbors
	const vector<int> TotalNeigh = calc_TotalNeigh(VolIndex);
#ifdef RICH_MPI
	vector<Primitive> MPIcells;
	vector<vector<double> > MPItracer;
	vector<vector<int> > MPI_AMR_Send; // the indeces in the Nghostpoints that I want to recv hydro from other procs
	CreateGetPrimitiveList(ToRemove,_tessellation.GetGhostIndeces(),n,MPI_AMR_Send);
	vector<int> ToRemoveReduced;
	if(n<static_cast<int>(ToRemove.size()))
	{
	  ToRemoveReduced.resize(static_cast<int>(ToRemove.size())-n);
		copy(ToRemove.begin()+n,ToRemove.end(),ToRemoveReduced.begin());
		//for(int i=0;i<static_cast<int>(ToRemoveReduced.size());++i)
		//	ToRemoveReduced[i]-=n;
	}
	GetAMRExtensive(MPIcells,MPItracer,_cells,tracer_,traceractive,MPI_AMR_Send,
		_tessellation.GetDuplicatedProcs(),_eos,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetGhostIndeces(),ToRemoveReduced);
#endif
	const vector<Conserved> c_temp = calc_c_temp(TotalNeigh,
						     VolIndex,
						     dv,
						     _cells,
						     #ifdef RICH_MPI
						     MPIcells,
						     #endif
						     ToRemove,
						     n);
	const vector<vector<double> > t_temp = calc_t_temp(traceractive,
							   TotalNeigh,
							   VolIndex,
							   dv,
							   _cells,
							   tracer_,
#ifdef RICH_MPI
							   MPIcells,
							   MPItracer,
#endif
							   ToRemove,
							   n);
	//	int Nneigh=int(TotalNeigh.size());
	// Update the primitives
	sort(ToRemove.begin(),ToRemove.end());
	remove_update_cells_tracers(traceractive,
				    TotalNeigh,
				    ToRemove,
				    _tessellation,
				    OldVol,
				    c_temp,
				    t_temp,
				    _eos,
				    _cells,
				    tracer_);
	if(n<static_cast<int>(ToRemove.size()))
		ToRemove.erase(ToRemove.begin()+n,ToRemove.end());
	//Remove the deleted cells
	RemoveVector(_cells,ToRemove);
	RemoveVector(custom_evolution_indices,ToRemove);
	//	RemoveVector(_conservedextensive,ToRemove);
	_conservedextensive = CalcConservedExtensive
	  (CalcConservedIntensive(_cells),_tessellation,*pg_);
	if(traceractive)
		RemoveVector(tracer_,ToRemove);
	// Fix the ghost points
	int nprocs=static_cast<int>(_tessellation.GetDuplicatedProcs().size());
	//sort(ToRemoveReduced.begin(),ToRemoveReduced.end());
	vector<vector<int> > & ghostpoints=_tessellation.GetDuplicatedPoints();
	#ifdef RICH_MPI
	vector<vector<int> > & nghost=_tessellation.GetGhostIndeces();
	vector<vector<int> > toremoveall(nprocs); // the indeces in the ghost that are removed
#endif // RICH_MPI
	for(int i=0;i<nprocs;++i)
	{
		vector<int> toremove2;
		int nsent2=static_cast<int>(ghostpoints[static_cast<size_t>(i)].size());
		for(int j=0;j<nsent2;++j)
		{
			int toReduce2=int(lower_bound(ToRemove.begin(),ToRemove.end(),ghostpoints[static_cast<size_t>(i)][static_cast<size_t>(j)])-
				ToRemove.begin());
			if(binary_search(ToRemove.begin(),ToRemove.end(),ghostpoints[static_cast<size_t>(i)][static_cast<size_t>(j)]))
				toremove2.push_back(j);
			else
				ghostpoints[static_cast<size_t>(i)][static_cast<size_t>(j)]-=toReduce2;
		}
#ifdef RICH_MPI
		if(!toremove2.empty())
			RemoveVector(ghostpoints[static_cast<size_t>(i)],toremove2);
		toremoveall[static_cast<size_t>(i)]=toremove2;
		nsent2=static_cast<int>(nghost[static_cast<size_t>(i)].size());
		for(int j=0;j<nsent2;++j)
			nghost[static_cast<size_t>(i)][static_cast<size_t>(j)]-=static_cast<int>(ToRemove.size());
#endif
	}
#ifdef RICH_MPI
	RemoveNGhostAMR(nghost,_tessellation.GetDuplicatedProcs(),toremoveall);
#endif
	// If self gravity is needed change here
	return ToRemove;
}

vector<int> hdsim::RefineCells(RefineStrategy *refine,vector<int>
	const& Removed,double dr)
{
  if(!refine)
    throw UniversalError("Error in refine, NULL pointer");
	const bool traceractive = !tracer_.empty();
	vector<Vector2D> directions;
	vector<int> PointsToRefine = refine->CellsToRefine(_tessellation,_cells,
		tracer_,_time,directions,Removed);
	// Get the list of points to refine
	if(PointsToRefine.empty())
		return PointsToRefine;
	PointsToRefine=refine->RemoveNearBoundary(PointsToRefine,directions,_tessellation);
	int n=int(PointsToRefine.size());
	const int N=_tessellation.GetPointNo();
	// Resize vectors
	_cells.resize(static_cast<size_t>(N+n));
	custom_evolution_indices.resize(static_cast<size_t>(N+n));
	if(traceractive)
	  tracer_.resize(static_cast<size_t>(N+n));
	// Change the mesh
	_tessellation.RefineCells(PointsToRefine,directions,dr);
	// Fill the new hydro
	for(int i=N;i<N+n;++i)
	{
	  _cells[static_cast<size_t>(i)]=_cells[static_cast<size_t>(PointsToRefine[static_cast<size_t>(i-N)])];
		if(traceractive)
		  tracer_[static_cast<size_t>(i)]=tracer_[static_cast<size_t>(PointsToRefine[static_cast<size_t>(i-N)])];
	}
	_conservedextensive = CalcConservedExtensive
	  (CalcConservedIntensive(_cells),_tessellation,*pg_);
	return PointsToRefine;
}

double hdsim::GetCfl(void)const
{
	return _cfl;
}

bool hdsim::GetColdFlowFlag(void)const
{
	return coldflows_flag_;
}

void hdsim::GetColdFlowParm(double &a,double &b)const
{
  a=cfp_.as;
  b=cfp_.bs;
}

bool hdsim::GetDensityFloorFlag(void) const
{
	return densityfloor_;
}

void hdsim::GetDensityFloorParm(double &density,double &pressure) const
{
	density=densityMin_;
	pressure=pressureMin_;
}

void hdsim::SetDensityFloor(double density,double pressure)
{
	densityfloor_=true;
	densityMin_=density;
	pressureMin_=pressure;
}

vector<Primitive>& hdsim::GetAllCells(void)
{
	return _cells;
}

vector<Conserved>& hdsim::getAllConserved(void)
{
  return _conservedextensive;
}

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

#ifdef RICH_MPI
void hdsim::SetProcessorMovement(ProcessorUpdate *procupdate)
{
	procupdate_=procupdate;
}
#endif

void hdsim::load(const ResetDump& checkpoint)
{
#ifdef RICH_MPI
	_proctess.Update(checkpoint.procmesh);
	_tessellation.Update(checkpoint.snapshot.mesh_points,_proctess);
#else
	_tessellation.Update(checkpoint.snapshot.mesh_points);
#endif
	_cells = checkpoint.snapshot.cells;
	_conservedextensive = CalcConservedExtensive
	  (CalcConservedIntensive(_cells),_tessellation,*pg_);
	_cfl = checkpoint.cfl;
	_time = checkpoint.time;
	cycle_ = checkpoint.cycle;
	tracer_ = checkpoint.tracers;
	tracer_flag_ = !tracer_.empty();
	coldflows_flag_=checkpoint.coldflows;
	densityfloor_ = checkpoint.densityfloor;
	cfp_.as = checkpoint.a;
	cfp_.bs = checkpoint.b;
	densityMin_ = checkpoint.densitymin;
	pressureMin_ = checkpoint.pressuremin;
}

void hdsim::makeCheckpoint(ResetDump& checkpoint) const
{
	checkpoint.snapshot.mesh_points = _tessellation.GetMeshPoints();
	int n=_tessellation.GetPointNo();
	checkpoint.snapshot.mesh_points.resize(static_cast<size_t>(n));
#ifdef RICH_MPI
	checkpoint.procmesh=_proctess.GetMeshPoints();
	n=_proctess.GetPointNo();
#endif
	checkpoint.procmesh.resize(static_cast<size_t>(n));
	checkpoint.snapshot.cells = _cells;
	checkpoint.cfl = _cfl;
	checkpoint.time = _time;
	checkpoint.cycle = cycle_;
	checkpoint.tracers = tracer_;
	checkpoint.coldflows = coldflows_flag_;
	checkpoint.densityfloor =  densityfloor_;
	checkpoint.a = cfp_.as;
	checkpoint.b = cfp_.bs;
	checkpoint.densitymin = densityMin_;
	checkpoint.pressuremin = pressureMin_;
}

void hdsim::setStartTime(double t_start)
{
	_time = t_start;
}

const EquationOfState& hdsim::getEos(void) const
{
  return _eos;
}

const OuterBoundary& hdsim::getOuterBoundary(void) const
{
  return _obc;
}

const HydroBoundaryConditions& hdsim::getHydroBoundaryCondition(void) const
{
  return _hbc;
}
