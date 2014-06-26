#include <cmath>
#include <algorithm>
#include "hdsim2d.hpp"
#include "hydrodynamics_2d.hpp"

using namespace std;

void hdsim::SetData(vector<Primitive> const& cells,
		    vector<Vector2D> const& points,
		    double time,
		    vector<vector<double> > const& tracers)
{
  _cells=cells;
  _tessellation.Update(points);
  tracer_=tracers;
  _time=time;
}

hdsim::hdsim
(vector<Vector2D> const& points,
 Tessellation& tessellation,
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
  _proctess(tessellation),
  _cells(vector<Primitive>()),
  _fluxes(vector<Conserved>()),
  _pointvelocity(vector<Vector2D>(points.size(),Vector2D(0,0))),
  _facevelocity(vector<Vector2D>()),
  _conservedintensive(vector<Conserved>()),
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
  as_(0),
  bs_(0),
  densityMin_(0),
  pressureMin_(0),
  EntropyReCalc_(EntropyCalc),
  _dt_external(0),
  custom_evolution_manager(),
  custom_evolution_indices(points.size(),0)
{
  _tessellation.Initialise(points, &obc);

  _cells = InitialiseCells(density, pressure,xvelocity, yvelocity,
			   eos, tessellation,CMvalue);

  _conservedintensive = CalcConservedIntensive(_cells);
  _conservedextensive = CalcConservedExtensive
    (_conservedintensive,tessellation);
  _dt_external=-1;
}

hdsim::hdsim
(vector<Vector2D> const& points,
 Tessellation& tessellation,
 Tessellation& proctess,
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
  _proctess(proctess),
  _cells(vector<Primitive>()),
  _fluxes(vector<Conserved>()),
  _pointvelocity(vector<Vector2D>(points.size(),Vector2D(0,0))),
  _facevelocity(vector<Vector2D>()),
  _conservedintensive(vector<Conserved>()),
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
  as_(0),
  bs_(0),
  densityMin_(0),
  pressureMin_(0),
  EntropyReCalc_(EntropyCalc),
  _dt_external(0),
  custom_evolution_manager(),
  custom_evolution_indices(points.size(),0)
{

#ifdef RICH_MPI
  int ws;
  MPI_Comm_size(MPI_COMM_WORLD,&ws);
  if(ws%2==1)
    throw UniversalError("MPI needs even number of threads");
#endif

  _tessellation.Initialise(points,_proctess,&obc);

  _cells = InitialiseCells(density, pressure,xvelocity, yvelocity,
			   eos, tessellation,CMvalue);

  _conservedintensive = CalcConservedIntensive(_cells);
  _conservedextensive = CalcConservedExtensive
    (_conservedintensive,tessellation);
  _dt_external=-1;
}

#ifdef RICH_MPI
hdsim::hdsim(ResetDump const& dump,Tessellation& tessellation,
	     Tessellation &tproc,
	     SpatialReconstruction& interpolation,
	     EquationOfState const& eos,RiemannSolver const& rs,
	     PointMotion& pointmotion,SourceTerm& external_force,
	     OuterBoundary const& obc,HydroBoundaryConditions const& hbc,
	     bool EntropyCalc):
  _tessellation(tessellation),
  _proctess(tproc),
  _cells(dump.snapshot.cells),
  _fluxes(vector<Conserved>()),
  _pointvelocity(vector<Vector2D>()),
  _facevelocity(vector<Vector2D>()),
  _conservedintensive(CalcConservedIntensive(_cells)),
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
  as_(dump.a),
  bs_(dump.b),
  densityMin_(dump.densitymin),
  pressureMin_(dump.pressuremin),
  EntropyReCalc_(EntropyCalc),
  _dt_external(-1),
  custom_evolution_manager(),
  custom_evolution_indices(dump.cevolve)
{

  int ws;
  MPI_Comm_size(MPI_COMM_WORLD,&ws);
  if(ws%2==1)
    throw UniversalError("MPI needs even number of threads");

  _tessellation.Initialise(dump.snapshot.mesh_points,_proctess,&_obc);
  _conservedextensive = CalcConservedExtensive
    (_conservedintensive,tessellation);
}
#endif


hdsim::hdsim(ResetDump const& dump,Tessellation& tessellation,
	     SpatialReconstruction& interpolation,
	     EquationOfState const& eos,RiemannSolver const& rs,
	     PointMotion& pointmotion,SourceTerm& external_force,
	     OuterBoundary const& obc,HydroBoundaryConditions const& hbc,
	     bool EntropyCalc):
  _tessellation(tessellation),
  _proctess(tessellation),
  _cells(dump.snapshot.cells),
  _fluxes(vector<Conserved>()),
  _pointvelocity(vector<Vector2D>()),
  _facevelocity(vector<Vector2D>()),
  _conservedintensive(CalcConservedIntensive(_cells)),
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
  as_(dump.a),
  bs_(dump.b),
  densityMin_(dump.densitymin),
  pressureMin_(dump.pressuremin),
  EntropyReCalc_(EntropyCalc),
  _dt_external(-1),
  custom_evolution_manager(),
  custom_evolution_indices(dump.cevolve)
{
  _tessellation.Initialise(dump.snapshot.mesh_points,&_obc);
  _conservedextensive = CalcConservedExtensive
    (_conservedintensive,tessellation);
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
}

void hdsim::TimeAdvance(void)
{

#ifndef RICH_MPI
  PeriodicUpdateCells(_cells,tracer_,custom_evolution_indices,
		      _tessellation.GetDuplicatedPoints(),_tessellation.GetTotalPointNumber());
#else
  SendRecvHydro(_cells,tracer_,custom_evolution_indices,_tessellation.GetDuplicatedPoints(),
		_tessellation.GetDuplicatedProcs(),_eos,_tessellation.GetTotalPointNumber());

  SendRecvVelocity(_tessellation.GetAllCM(),_tessellation.GetDuplicatedPoints(),
		   _tessellation.GetDuplicatedProcs(),_tessellation.GetTotalPointNumber());
#endif

  vector<CustomEvolution*> custom_evolutions = 
    convert_indices_to_custom_evolution(custom_evolution_manager,
					custom_evolution_indices);

  CalcPointVelocities(_tessellation, _cells, 
		      _pointmotion, _pointvelocity,_time,custom_evolutions);

#ifndef RICH_MPI
  PeriodicVelocityExchange(_pointvelocity,_tessellation.GetDuplicatedPoints(),
			   _tessellation.GetTotalPointNumber());
#else
  SendRecvVelocity(_pointvelocity,_tessellation.GetDuplicatedPoints(),
		   _tessellation.GetDuplicatedProcs(),_tessellation.GetTotalPointNumber());
#endif

  _facevelocity=_tessellation.calc_edge_velocities
    (&_hbc,_pointvelocity,_time);

  double dt = _cfl*CalcTimeStep(_tessellation, _cells, _facevelocity,_hbc,
				_time,custom_evolutions);

  if(_dt_external>0)
    dt=min(dt,_dt_external*_cfl);

  if(_endtime>0)
    if(_time+dt>_endtime)
      dt=_endtime-_time;

#ifdef RICH_MPI
  double dt_temp=dt;
  MPI_Reduce(&dt_temp,&dt,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

  vector<double> Ek,Ef;
  vector<char> shockedcells;
  if(coldflows_flag_)
    {
      const int n=_tessellation.GetPointNo();
      for(int i=0;i<n;++i)
	tracer_[i][0]=_eos.dp2s(_cells[i].Density,_cells[i].Pressure);
      shockedcells.resize(n);
      for(int i=0;i<n;++i)
	shockedcells[i]=IsShockedCell(_tessellation,i,_cells,_hbc,_time) ? 1 : 0 ;
    }

  CalcFluxes(_tessellation, _cells, dt, _time,
	     _interpolation,
	     _facevelocity, _hbc, _rs,
	     _fluxes,custom_evolutions,
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
	(tracer_,_cells,_tessellation,_fluxes,dt,_hbc,
	 _interpolation,_time,custom_evolutions,
	 custom_evolution_manager,
	 _facevelocity,lengths);
      MakeTracerExtensive(tracer_,
			  _tessellation,_cells,tracer_extensive);
      UpdateTracerExtensive(tracer_extensive,
			    trace_change,custom_evolutions,
			    _cells, _tessellation,_time);
    }

  _conservedintensive=CalcConservedIntensive(_cells);

  _conservedextensive=CalcConservedExtensive(_conservedintensive,_tessellation);

  vector<double> g;
  ExternalForceContribution(_tessellation,_cells,external_force_,_time,dt,
			    _conservedextensive,_hbc,_fluxes,_pointvelocity,g,coldflows_flag_,tracer_,
			    lengths);

  if(coldflows_flag_)
    {
      Ek =GetMaxKineticEnergy(_tessellation,_cells,custom_evolutions);
      Ef = GetForceEnergy(_tessellation,g);
    }


  UpdateConservedExtensive(_tessellation, _fluxes, dt,
			   _conservedextensive,_hbc,lengths);

#ifndef RICH_MPI
  MoveMeshPoints(_pointvelocity, dt, _tessellation);
#else
  if(procupdate_!=0)
    procupdate_->Update(_proctess,_tessellation);
  MoveMeshPoints(_pointvelocity, dt, _tessellation,_proctess);
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

  if(!ptoadd.empty())
    {
      _conservedextensive.insert(_conservedextensive.end(),ptoadd.begin(),
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


  custom_evolutions = 
    convert_indices_to_custom_evolution(custom_evolution_manager,
					custom_evolution_indices);

  vector<char> btoadd;
  if(coldflows_flag_)
    {
      vector<char> btoadd;
      SendRecvShockedStatus(shockedcells,_tessellation.GetSentPoints(),
			    _tessellation.GetSentProcs(),btoadd);
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

  UpdateConservedIntensive(_tessellation, _conservedextensive, 
			   _conservedintensive);

  if(coldflows_flag_)
    {
      FixPressure(_conservedintensive,tracer_extensive,_eos,Ek,Ef,as_,bs_,
		  custom_evolutions,_tessellation,_conservedextensive,shockedcells,
		  densityfloor_);
    }

  UpdatePrimitives(_conservedintensive, _eos, _cells,custom_evolutions,_cells,
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

Tessellation const& hdsim::GetProcTessellation(void) const
{
  return _proctess;
}

void hdsim::TimeAdvance2Mid(void)
{
  const double dt=TimeAdvance2mid
    (_tessellation,_proctess,_cells,_pointmotion,
     _hbc,_interpolation,_rs,_eos,external_force_,_time,_cfl,_endtime,
     tracer_,_dt_external,custom_evolution_indices,
     custom_evolution_manager,
     #ifdef RICH_MPI
     procupdate_,
     #endif
     tracer_flag_,
     coldflows_flag_,as_,bs_,densityfloor_,densityMin_,pressureMin_,
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
    tracer_[i].push_back(tp.EvalAt(_tessellation.GetCellCM(i)));
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

Conserved hdsim::GetFlux(int i) const
{
  return _fluxes[i];
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
  as_=as;
  bs_=bs;
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
			SpatialDistribution const& originalVy,int tracerindex)
{
  vector<CustomEvolution*> cevolve=convert_indices_to_custom_evolution(
								       custom_evolution_manager,custom_evolution_indices);
  TracerResetCalc(alpha,originalD,originalP,originalVx,originalVy,
		  _cells,_tessellation,tracer_,tracerindex,_eos,cevolve);
  return;
}


vector<int> hdsim::RemoveCells(RemovalStrategy const* remove)
{
  if(!remove)
    throw UniversalError("No Removal strategy");
  vector<int> ToRemove=remove->CellsToRemove(_tessellation,_cells,tracer_,_time);
  if(ToRemove.empty())
    return ToRemove;
  int n=int(ToRemove.size());
  bool traceractive;
  if(!tracer_.empty())
    traceractive=true;
  else
    traceractive=false;
  sort(ToRemove.begin(),ToRemove.end());
  // save the extensive of the removed cells
  vector<double> OldVol(_tessellation.GetPointNo());
  vector<Conserved> OldExtensive;
  vector<vector<double> > OldTracer;
  OldExtensive.reserve(n);
  if(traceractive)
    OldTracer.reserve(n);
  for(int i=0;i<(int)_tessellation.GetPointNo();++i)
    OldVol[i]=_tessellation.GetVolume(i);
  for(int i=0;i<n;++i)
    {
      double vol=OldVol[ToRemove[i]];
      OldExtensive.push_back(Primitive2Conserved(_cells[ToRemove[i]],vol));
      if(traceractive)
	OldTracer.push_back(_cells[ToRemove[i]].Density*vol*
			    tracer_[ToRemove[i]]);
    }
  // Change the tessellation
  vector<vector<int> > VolIndex;
  vector<vector<double> > VolRatio;
  _tessellation.RemoveCells(ToRemove,VolIndex,VolRatio);
  // gather all the relevant neighbors
  vector<int> TotalNeigh;
  for(int i=0;i<n;++i)
    TotalNeigh.insert(TotalNeigh.end(),VolIndex[i].begin(),
		      VolIndex[i].end());
  sort(TotalNeigh.begin(),TotalNeigh.end());
  TotalNeigh=unique(TotalNeigh);
  // Calculate the old extensive hydro
  vector<Conserved> c_temp;
  vector<vector<double> > t_temp;
  int Nneigh=int(TotalNeigh.size());
  c_temp.reserve(Nneigh);
  if(traceractive)
    t_temp.reserve(Nneigh);
  for(int i=0;i<Nneigh;++i)
    {
      c_temp.push_back(Primitive2Conserved(_cells[TotalNeigh[i]],
					   OldVol[TotalNeigh[i]]));
      if(traceractive)
	t_temp.push_back(_cells[TotalNeigh[i]].Density*
			 OldVol[TotalNeigh[i]]*tracer_[TotalNeigh[i]]);
    }
  // Change the extensive
  for(int i=0;i<n;++i)
    {
      for(int j=0;j<(int)VolIndex[i].size();++j)
	{
	  int index=int(lower_bound
			(TotalNeigh.begin(),TotalNeigh.end(),
			 VolIndex[i][j])-TotalNeigh.begin());
	  c_temp[index]+=VolRatio[i][j]*OldExtensive[i];
	  if(traceractive)
	    for(int jj=0;jj<(int)t_temp[0].size();++jj)
	      t_temp[index][jj]+=VolRatio[i][j]*OldTracer[i][jj];
	}
    }
  // Update the primitives
  for(int i=0;i<Nneigh;++i)
    {
      int index=int(lower_bound(ToRemove.begin(),ToRemove.end(),TotalNeigh[i])-
		    ToRemove.begin());
      double vol_inv=1.0/_tessellation.GetVolume(TotalNeigh[i]-index);
      try
	{
	  _cells[TotalNeigh[i]]=Conserved2Primitive(vol_inv*c_temp[i],_eos);
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
	  tracer_[TotalNeigh[i]]=density_inv*vol_inv*t_temp[i];
	}
    }
  //Remove the deleted cells
  RemoveVector(_cells,ToRemove);
  RemoveVector(custom_evolution_indices,ToRemove);
  if(traceractive)
    RemoveVector(tracer_,ToRemove);
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
  a=as_;
  b=bs_;
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

vector<Conserved>const& hdsim::GetFluxes(void)const
{
  return _fluxes;
}

Vector2D hdsim::GetPointVelocity(int index)const
{
  return _pointvelocity[index];
}

vector<Vector2D> const& hdsim::getAllPointVelocities(void) const
{
  return _pointvelocity;
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
  _conservedintensive = CalcConservedIntensive(_cells);
  _conservedextensive = CalcConservedExtensive
    (_conservedintensive,_tessellation);
  _cfl = checkpoint.cfl;
  _time = checkpoint.time;
  cycle_ = checkpoint.cycle;
  tracer_ = checkpoint.tracers;
  tracer_flag_ = !tracer_.empty();
  coldflows_flag_=checkpoint.coldflows;
  densityfloor_ = checkpoint.densityfloor;
  as_ = checkpoint.a;
  bs_ = checkpoint.b;
  densityMin_ = checkpoint.densitymin;
  pressureMin_ = checkpoint.pressuremin;
}

void hdsim::makeCheckpoint(ResetDump& checkpoint) const
{
  checkpoint.snapshot.mesh_points = _tessellation.GetMeshPoints();
  int n=_tessellation.GetPointNo();
  checkpoint.snapshot.mesh_points.resize(n);
  checkpoint.procmesh=_proctess.GetMeshPoints();
  n=_proctess.GetPointNo();
  checkpoint.procmesh.resize(n);
  checkpoint.snapshot.cells = _cells;
  checkpoint.cfl = _cfl;
  checkpoint.time = _time;
  checkpoint.cycle = cycle_;
  checkpoint.tracers = tracer_;
  checkpoint.coldflows = coldflows_flag_;
  checkpoint.densityfloor =  densityfloor_;
  checkpoint.a = as_;
  checkpoint.b = bs_;
  checkpoint.densitymin = densityMin_;
  checkpoint.pressuremin = pressureMin_;
}

void hdsim::setStartTime(double t_start)
{
  _time = t_start;
}
