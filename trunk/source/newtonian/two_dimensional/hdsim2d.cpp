#include <cmath>
#include <algorithm>
#include "hdsim2d.hpp"
#include "hydrodynamics_2d.hpp"

using namespace std;

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
	  (CalcConservedIntensive(_cells),tessellation);
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
	  (CalcConservedIntensive(_cells),tessellation);
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
    vector<char> res(tess.GetPointNo());
    for(size_t i=0;i<res.size();++i)
      res[i] = IsShockedCell(tess,i,cells,hbc,time) ? 1 : 0;
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

    const vector<char>& getShockedStatus(void)
    {
      return shocked_cells_;
    }

    const vector<double>& getKineticEnergy(void)
    {
      return kinetic_energy_;
    }

    const vector<double>& getPotentialEnergy(void)
    {
      return potential_energy_;
    }

    void fixPressure(vector<Conserved>& intensive,
		     vector<vector<double> >& tracers,
		     const EquationOfState& eos,
		     double as,
		     double bs,
		     vector<CustomEvolution*>& ce,
		     Tessellation& tess,
		     vector<Conserved>& extensive,
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

  void substitute_entropy(const vector<Primitive>& cells,
			  const EquationOfState& eos,
			  size_t index,
			  vector<vector<double> >& tracers)
  {
    for(size_t i=0;i<cells.size();++i)
      tracers.at(i).at(index) = eos.dp2s(cells.at(i).Density,cells.at(i).Pressure);
  }
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
	  substitute_entropy(_cells,_eos,0,tracer_);

#ifdef RICH_MPI
	if(tracer_flag_)
		SendRecvTracers(tracer_,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_tessellation.GetGhostIndeces()
		,_tessellation.GetTotalPointNumber());
#endif
	vector<Conserved> fluxes = calc_fluxes
	  (_tessellation, _cells, dt, _time,
	   _interpolation,
	   fv, _hbc, _rs,
	   custom_evolutions,
	   custom_evolution_manager,
	   tracer_);
						

	int nn=_tessellation.GetTotalSidesNumber();
	vector<double> lengths(nn);
	for(int i=0;i<nn;++i)
		lengths[i]=_tessellation.GetEdge(i).GetLength();

	vector<vector<double> > tracer_extensive;
	if(tracer_flag_)
	{
		vector<vector<double> > trace_change;
		trace_change = CalcTraceChange
			(tracer_,_cells,_tessellation,fluxes,dt,_hbc,
			_interpolation,_time,custom_evolutions,
			custom_evolution_manager,
			fv,lengths);
		MakeTracerExtensive(tracer_,
			_tessellation,_cells,tracer_extensive);
		UpdateTracerExtensive(tracer_extensive,
			trace_change,custom_evolutions,
			_cells, _tessellation,_time);
	}

	/*
	_conservedintensive=CalcConservedIntensive(_cells);

	_conservedextensive=CalcConservedExtensive(_conservedintensive,_tessellation);
	*/

	vector<double> g;
	ExternalForceContribution(_tessellation,_cells,external_force_,_time,dt,
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
	if(procupdate_!=0)
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

	vector<Conserved> intensive;
	UpdateConservedIntensive
	  (_tessellation, _conservedextensive, intensive);

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

	if(tracer_flag_)
	{
		MakeTracerIntensive(tracer_,tracer_extensive,_tessellation,_cells);
	}

	_time += dt;
	cycle_++;
}

Tessellation const& hdsim::GetTessellation(void) const
{
	return _tessellation;
}

#ifdef RICH_MPI
Tessellation const& hdsim::GetProcTessellation(void) const
{
	return _proctess;
}
#endif

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
		custom_evolution_manager,
#ifdef RICH_MPI
		procupdate_,
#endif
		tracer_flag_,
		coldflows_flag_,cfp_.as,cfp_.bs,densityfloor_,densityMin_,pressureMin_,
		EntropyReCalc_);

	_time += dt;
	++cycle_;
}

void hdsim::addTracer(SpatialDistribution const& tp)
{
	const int n = _tessellation.GetPointNo();
	if(!tracer_flag_){
		tracer_.resize(n);
		tracer_flag_ = true;
	}

	for(int i=0;i<n;++i)
		tracer_[i].push_back(tp(_tessellation.GetCellCM(i)));
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
	return _cells[i];
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
		tracer_.resize(n);
		tracer_flag_ = true;
	}

	for(int i=0;i<n;++i)
		tracer_[i].push_back(0);
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
	void CreateGetPrimitiveList(vector<int> const& ToRemove,vector<vector<int> >
		const& Nghost,int nremoved,vector<vector<int> > &MPI_AMR_Send)
	{
		int nprocs=(int)Nghost.size();
		MPI_AMR_Send.resize(nprocs);
		vector<vector<int> > SortedNghost(nprocs),SortIndex(nprocs);
		// sort Nghost
		for(int i=0;i<nprocs;++i)
		{
			SortedNghost[i]=Nghost[i];
			sort_index(Nghost[i],SortIndex[i]);
			sort(SortedNghost[i].begin(),SortedNghost[i].end());
		}
		for(int i=0;i<(int)ToRemove.size();++i)
		{
			for(int j=0;j<nprocs;++j)
			{
				if(binary_search(SortedNghost[j].begin(),SortedNghost[j].end(),
					ToRemove[i]))
				{
					int index2=lower_bound(SortedNghost[j].begin(),SortedNghost[j].end(),
						ToRemove[i])-SortedNghost[j].begin();
					MPI_AMR_Send[j].push_back(SortIndex[j][index2]);
				}
			}
		}
	}
}

namespace
{
#ifdef RICH_MPI
	void RemoveNGhostAMR(vector<vector<int> > &nghost,vector<int> const& sentprocs,
		vector<vector<int> > &toremove)
	{
		int nlist=(int)sentprocs.size();
		int rank=get_mpi_rank();
		int ws=get_mpi_size();
		vector<int> procorder=GetProcOrder(rank,ws);
		vector<vector<int> > recv(nlist);
		int temp;
		MPI_Status status;
		for(int i=0;i<(int)procorder.size();++i)
		{
			int index=Find(sentprocs.begin(),sentprocs.end(),procorder[i])
				-sentprocs.begin();
			if(index<nlist)
			{
				if(rank<procorder[i])
				{
					if(toremove[index].empty())
						MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
					else
						MPI_Send(&toremove[index][0],(int)toremove[index].size(),
						MPI_INT,procorder[i],0,MPI_COMM_WORLD);
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
						MPI_Send(&toremove[index][0],(int)toremove[index].size(),
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

vector<int> hdsim::RemoveCells(RemovalStrategy const* remove)
{
	if(!remove)
		throw UniversalError("No Removal strategy");
	vector<int> ToRemove=remove->CellsToRemove(_tessellation,_cells,tracer_,_time);
	//if(ToRemove.empty())
	//	return ToRemove;
	int n=int(ToRemove.size());
	bool traceractive;
	if(!tracer_.empty())
		traceractive=true;
	else
		traceractive=false;
	if(!ToRemove.empty())
		sort(ToRemove.begin(),ToRemove.end());
	// save the extensive of the removed cells
	vector<double> OldVol(_tessellation.GetPointNo());
	for(int i=0;i<(int)_tessellation.GetPointNo();++i)
		OldVol[i]=_tessellation.GetVolume(i);
	// Change the tessellation
	vector<vector<int> > VolIndex;
	vector<vector<double> > dv;
	_tessellation.RemoveCells(ToRemove,VolIndex,dv);
	n=lower_bound(ToRemove.begin(),ToRemove.end(),(int)OldVol.size())-
		ToRemove.begin();
	// gather all the relevant neighbors
	vector<Primitive> MPIcells;
	vector<vector<double> > MPItracer;
	vector<int> TotalNeigh;
	for(int i=0;i<(int)dv.size();++i)
		TotalNeigh.insert(TotalNeigh.end(),VolIndex[i].begin(),
		VolIndex[i].end());
	sort(TotalNeigh.begin(),TotalNeigh.end());
	TotalNeigh=unique(TotalNeigh);
#ifdef RICH_MPI
	int rank=get_mpi_rank();
	vector<vector<int> > MPI_AMR_Send; // the indeces in the Nghostpoints that I want to recv hydro from other procs
	CreateGetPrimitiveList(ToRemove,_tessellation.GetGhostIndeces(),n,MPI_AMR_Send);
	vector<int> ToRemoveReduced;
	if(n<(int)ToRemove.size())
	{
		ToRemoveReduced.resize((int)ToRemove.size()-n);
		copy(ToRemove.begin()+n,ToRemove.end(),ToRemoveReduced.begin());
		//for(int i=0;i<(int)ToRemoveReduced.size();++i)
		//	ToRemoveReduced[i]-=n;
	}
	GetAMRExtensive(MPIcells,MPItracer,_cells,tracer_,traceractive,MPI_AMR_Send,
		_tessellation.GetDuplicatedProcs(),_eos,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetGhostIndeces(),ToRemoveReduced);
#endif
	// Calculate the old extensive hydro
	vector<Conserved> c_temp;
	vector<vector<double> > t_temp;
	int Nneigh=int(TotalNeigh.size());
	c_temp.resize(Nneigh); // need to add hydro/tracer of removed ghost
	if(traceractive)
	{
		t_temp.resize(Nneigh);
		for(int i=0;i<Nneigh;++i)
		{
			t_temp[i].resize(tracer_[0].size());
		}
	}
	// Change the extensive
	for(int i=0;i<(int)dv.size();++i)
	{
		for(int j=0;j<(int)VolIndex[i].size();++j)
		{
			int index=int(lower_bound(TotalNeigh.begin(),TotalNeigh.end(),
				VolIndex[i][j])-TotalNeigh.begin());
			if(i<n)
				c_temp[index]+=Primitive2Conserved(_cells[ToRemove[i]],dv[i][j]);
			else
			{
				c_temp[index]+=Primitive2Conserved(MPIcells[i-n],dv[i][j]);
			}
			if(traceractive)
				if(i<n)
					for(int jj=0;jj<(int)t_temp[0].size();++jj)
						t_temp[index][jj]+=dv[i][j]*tracer_[i][jj]*_cells[ToRemove[i]].Density;
				else
					for(int jj=0;jj<(int)t_temp[0].size();++jj)
						t_temp[index][jj]+=dv[i][j]*MPItracer[i-n][jj]*MPIcells[i-n].Density;
		}
	}
	// Update the primitives
	sort(ToRemove.begin(),ToRemove.end());
	for(int i=0;i<Nneigh;++i)
	{
		int index=int(lower_bound(ToRemove.begin(),ToRemove.end(),TotalNeigh[i])-
			ToRemove.begin());
		double vol_inv=1.0/_tessellation.GetVolume(TotalNeigh[i]-index);
		Conserved oldextensive=Primitive2Conserved(_cells[TotalNeigh[i]],
			OldVol[TotalNeigh[i]]);
		double olddensity=_cells[TotalNeigh[i]].Density;
		try
		{
			_cells[TotalNeigh[i]]=Conserved2Primitive(vol_inv*(c_temp[i]+
				oldextensive),_eos);
		}
		catch(UniversalError &eo)
		{
			eo.AddEntry("Error in cell derefine",0);
			eo.AddEntry("volume of adjusted cell",1.0/vol_inv);
			eo.AddEntry("Index of adjusted cell",(double)TotalNeigh[i]);
			eo.AddEntry("Index of adjusted cell in Totalneigh",(double)i);
			eo.AddEntry("ToRemove size",(double)ToRemove.size());
			eo.AddEntry("Index of adjusted cell in ToRemove",(double)index);
			eo.AddEntry("Totalneigh size",(double)TotalNeigh.size());
			throw;
		}
		if(traceractive)
		{
			double density_inv=1.0/_cells[TotalNeigh[i]].Density;
			tracer_[TotalNeigh[i]]=density_inv*vol_inv*(t_temp[i]+
				OldVol[TotalNeigh[i]]*olddensity*tracer_[TotalNeigh[i]]);
		}
	}
	if(n<(int)ToRemove.size())
		ToRemove.erase(ToRemove.begin()+n,ToRemove.end());
	//Remove the deleted cells
	RemoveVector(_cells,ToRemove);
	RemoveVector(custom_evolution_indices,ToRemove);
	if(traceractive)
		RemoveVector(tracer_,ToRemove);
	// Fix the ghost points
	int nprocs=(int)_tessellation.GetDuplicatedProcs().size();
	//sort(ToRemoveReduced.begin(),ToRemoveReduced.end());
	vector<vector<int> > & ghostpoints=_tessellation.GetDuplicatedPoints();
	vector<vector<int> > & nghost=_tessellation.GetGhostIndeces();
	vector<vector<int> > toremoveall(nprocs); // the indeces in the ghost that are removed
	for(int i=0;i<nprocs;++i)
	{
		vector<int> toremove2;
		int nsent2=(int)ghostpoints[i].size();
		for(int j=0;j<nsent2;++j)
		{
			int toReduce2=int(lower_bound(ToRemove.begin(),ToRemove.end(),ghostpoints[i][j])-
				ToRemove.begin());
			if(binary_search(ToRemove.begin(),ToRemove.end(),ghostpoints[i][j]))
				toremove2.push_back(j);
			else
				ghostpoints[i][j]-=toReduce2;
		}
#ifdef RICH_MPI
		if(!toremove2.empty())
			RemoveVector(ghostpoints[i],toremove2);
		toremoveall[i]=toremove2;
		nsent2=(int)nghost[i].size();
		for(int j=0;j<nsent2;++j)
			nghost[i][j]-=(int)ToRemove.size();
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
	bool traceractive;
	if(!tracer_.empty())
		traceractive=true;
	else
		traceractive=false;
	vector<int> PointsToRefine;
	vector<Vector2D> directions;
	// Get the list of points to refine
	if(refine!=0)
		PointsToRefine=refine->CellsToRefine(_tessellation,_cells,
		tracer_,_time,directions,Removed);
	else
		throw UniversalError("Error in refine, NULL pointer");
	if(PointsToRefine.empty())
		return PointsToRefine;
	PointsToRefine=refine->RemoveNearBoundary(PointsToRefine,directions,_tessellation);
	int n=int(PointsToRefine.size());
	const int N=_tessellation.GetPointNo();
	// Resize vectors
	_cells.resize(N+n);
	custom_evolution_indices.resize(N+n);
	if(traceractive)
		tracer_.resize(N+n);
	// Change the mesh
	_tessellation.RefineCells(PointsToRefine,directions,dr);
	// Fill the new hydro
	for(int i=N;i<N+n;++i)
	{
		_cells[i]=_cells[PointsToRefine[i-N]];
		if(traceractive)
			tracer_[i]=tracer_[PointsToRefine[i-N]];
	}
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

void hdsim::HilbertArrange(int innernum)
{
	vector<Vector2D> cor=_tessellation.GetMeshPoints();
	vector<int> order=HilbertOrder(cor,_tessellation.GetPointNo(),innernum);
	ReArrangeVector(cor,order);
	if(cor.size()>order.size())
		cor.erase(cor.begin()+order.size(),cor.end());
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
	  (CalcConservedIntensive(_cells),_tessellation);
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
	checkpoint.snapshot.mesh_points.resize(n);
#ifdef RICH_MPI
	checkpoint.procmesh=_proctess.GetMeshPoints();
	n=_proctess.GetPointNo();
#endif
	checkpoint.procmesh.resize(n);
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
