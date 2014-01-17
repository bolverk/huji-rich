#include <cmath>
#include <algorithm>
#include <iostream>
#include "hdsim2d.hpp"
#include "hydrodynamics_2d.hpp"

using namespace std;

void hdsim::SetData(vector<Primitive> const& cells,
	vector<Vector2D> const& points,double time,vector<vector<double> > const& tracers)
{
	_cells=cells;
	_tessellation->Update(points);
	tracer_=tracers;
	_time=time;
}

hdsim::hdsim
(vector<Vector2D> const& points,
 Tessellation* tessellation,
 SpatialReconstruction* interpolation,
 SpatialDistribution const& density,
 SpatialDistribution const& pressure,
 SpatialDistribution const& xvelocity,
 SpatialDistribution const& yvelocity,
 EquationOfState const& eos,
 RiemannSolver const& rs,
 PointMotion *pointmotion,
 SourceTerm *external_force,
 OuterBoundary const* obc,
 HydroBoundaryConditions const* hbc,
 bool EntropyCalc,bool CMvalue):
  _tessellation(tessellation),
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
  CellsEvolve(points.size(),0)
{
  _tessellation->Initialise(points, obc);

  _cells = InitialiseCells(density, pressure,xvelocity, yvelocity,
			   eos, tessellation,CMvalue);

  _conservedintensive = CalcConservedIntensive(_cells);
  _conservedextensive = CalcConservedExtensive
    (_conservedintensive,tessellation);
  _dt_external=-1;
}

hdsim::hdsim(vector<Vector2D> const& points,Tessellation* tessellation,
	     SpatialReconstruction* interpolation,vector<Primitive> const& cells,
	     EquationOfState const& eos,RiemannSolver const& rs,
	     PointMotion *pointmotion,SourceTerm *external_force,
	     OuterBoundary const* obc,HydroBoundaryConditions const* hbc,
	     vector<vector<double> > const& tracers,double time,double cfl,
	     int cycle,bool coldflows,double a,double b,bool densityfloor,
	     double densitymin,double pressuremin,bool EntropyCalc):
  _tessellation(tessellation),
  _cells(cells),
  _fluxes(vector<Conserved>()),
  _pointvelocity(vector<Vector2D>()),
  _facevelocity(vector<Vector2D>()),
  _conservedintensive(vector<Conserved>()),
  _conservedextensive(vector<Conserved>()),
  _eos(eos),
  _rs(rs),
  _interpolation(interpolation),
  _pointmotion(pointmotion),
  _hbc(hbc),
  _obc(obc),
  external_force_(external_force),
  _cfl(cfl),
  _time(time),
  _endtime(-1),
  cycle_(cycle),
  tracer_(vector<vector<double> >()),
  tracer_flag_(false),
  coldflows_flag_(coldflows),
  densityfloor_(densityfloor),
  as_(a),
  bs_(b),
  densityMin_(densitymin),
  pressureMin_(pressuremin),
  EntropyReCalc_(EntropyCalc),
  _dt_external(0),
  CellsEvolve(vector<CustomEvolution*>())
{
	tracer_=tracers;			
	if(!tracer_.empty())
		tracer_flag_=true;
	else
	{
		tracer_flag_=false;
	}
	_tessellation->Initialise(points,_obc);
	_conservedintensive = CalcConservedIntensive(_cells);
	_conservedextensive = CalcConservedExtensive
		(_conservedintensive,tessellation);
	int n=int(_cells.size());
	CellsEvolve.reserve(n);
	for(int i=0;i<n;++i)
		CellsEvolve.push_back(0);
	_dt_external=-1;
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

void hdsim::TimeAdvance(void)
{

	CalcPointVelocities(_tessellation, _cells, 
		_pointmotion, _pointvelocity,_time);

	_facevelocity=_tessellation->calc_edge_velocities(_hbc,_pointvelocity,_time);
	
	double dt = _cfl*CalcTimeStep(_tessellation, _cells, _facevelocity,_hbc,
		_time);

	if(_dt_external>0)
		dt=min(dt,_dt_external*_cfl);

	if(_endtime>0)
		if(_time+dt>_endtime)
			dt=_endtime-_time;

	if(coldflows_flag_)
	{
	  int n=int(_cells.size());
		for(int i=0;i<n;++i)
			tracer_[i][0]=_eos.dp2s(_cells[i].Density,_cells[i].Pressure);
	}

	CalcFluxes(_tessellation, _cells, dt, _time,
		_interpolation,
		_facevelocity, _hbc, _rs,
		_fluxes,CellsEvolve,tracer_);

	vector<vector<double> > tracer_extensive;
	if(tracer_flag_)
	{
	  vector<vector<double> > trace_change;
	  trace_change = CalcTraceChange(tracer_,_cells,_tessellation,_fluxes,dt,_hbc,_interpolation,_time);
	  MakeTracerExtensive(tracer_,_tessellation,_cells,tracer_extensive);
	  UpdateTracerExtensive(tracer_extensive,trace_change);
	}

	_conservedintensive=CalcConservedIntensive(_cells);

	_conservedextensive=CalcConservedExtensive(_conservedintensive,_tessellation);

	UpdateConservedExtensive(_tessellation, _fluxes, dt,
		_conservedextensive,_hbc);

	vector<double> g;
	ExternalForceContribution
	  (_tessellation,_cells,
	   external_force_,
	   _time,
	   dt,
	   _conservedextensive,
	   _hbc,
	   _fluxes,
	   _pointvelocity,
	   g,
	   coldflows_flag_,tracer_);

	MoveMeshPoints(_pointvelocity, dt, _tessellation);

	UpdateConservedIntensive(_tessellation, _conservedextensive, 
		_conservedintensive);

	if(coldflows_flag_)
	{
	  const vector<double> Ek = 
	    GetMaxKineticEnergy(_tessellation,_cells,CellsEvolve);
	  const vector<double> Ef = GetForceEnergy(_tessellation,g);
	  FixPressure(_conservedintensive,tracer_extensive,_eos,Ek,Ef,as_,bs_,
		      CellsEvolve,_tessellation,_conservedextensive,_cells,_hbc,
		      _time);
	}

	UpdatePrimitives(_conservedintensive, _eos, _cells,CellsEvolve,_cells);

	if(tracer_flag_)
	{
		MakeTracerIntensive(tracer_,tracer_extensive,_tessellation,_cells);
	}

	_time += dt;
	cycle_++;	
}

Tessellation const* hdsim::GetTessellation(void) const
{
	return _tessellation;
}

void hdsim::TimeAdvance2Mid(void)
{
  const double dt=TimeAdvance2mid
    (_tessellation,_cells,_pointmotion,
     _hbc,_interpolation,_rs,_eos,external_force_,_time,_cfl,_endtime,
     CellsEvolve,tracer_,_dt_external,tracer_flag_,
     coldflows_flag_,as_,bs_,densityfloor_,densityMin_,pressureMin_,
     EntropyReCalc_);

  _time += dt;
  ++cycle_;
}

void hdsim::addTracer(SpatialDistribution const& tp)
{
	if(!tracer_flag_){
		tracer_.resize(_cells.size());
		tracer_flag_ = true;
	}

	for(int i=0;i<(int)_cells.size();++i)
		tracer_[i].push_back(tp.EvalAt(_tessellation->GetCellCM(i)));
}

// Diagnostics

int hdsim::GetEdgeNo(void) const
{
  return _tessellation->GetTotalSidesNumber();
}

Edge hdsim::GetEdge(int i) const
{
  return _tessellation->GetEdge(i);
}

int hdsim::GetCellNo(void) const
{
  return _tessellation->GetPointNo();
}

Primitive hdsim::GetCell(int i) const
{
  return _cells[i];
}

Vector2D hdsim::GetMeshPoint(int i) const
{
  return _tessellation->GetMeshPoint(i);
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
  return _tessellation->GetVolume(index);
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
	coldflows_flag_=true;
	if(!tracer_flag_)
	{
		tracer_.resize(_cells.size());
		tracer_flag_ = true;
	}

	for(int i=0;i<(int)_cells.size();++i)
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
	SpatialDistribution const& originalVy,int tracerindex,int InnerNumber)
{
	TracerResetCalc(alpha,originalD,originalP,originalVx,originalVy,
		_cells,_tessellation,tracer_,tracerindex,_eos,InnerNumber);
	return;
}


vector<int> hdsim::RemoveCells(RemovalStrategy const* remove)
{
	if(remove==0)
		throw UniversalError("No Removal strategy");
	vector<int> ToRemove=remove->CellsToRemove(_tessellation,_cells,tracer_,
		_time);
	int n=int(ToRemove.size());
	if(n==0)
		return ToRemove;
	bool traceractive;
	if(!tracer_.empty())
		traceractive=true;
	else
		traceractive=false;
	sort(ToRemove.begin(),ToRemove.end());
	// save the extensive of the removed cells
	vector<double> OldVol(_cells.size());
	vector<Conserved> OldExtensive;
	vector<vector<double> > OldTracer;
	OldExtensive.reserve(n);
	if(traceractive)
		OldTracer.reserve(n);
	for(int i=0;i<(int)_cells.size();++i)
		OldVol[i]=_tessellation->GetVolume(i);
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
	_tessellation->RemoveCells(ToRemove,VolIndex,VolRatio);
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
		double vol_inv=1.0/_tessellation->GetVolume(TotalNeigh[i]-index);
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
	RemoveVector(CellsEvolve,ToRemove);
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
	int n=int(PointsToRefine.size());
	int N=int(_cells.size());
	// Resize vectors
	_cells.resize(N+n);
	CellsEvolve.resize(N+n);
	if(traceractive)
		tracer_.resize(N+n);
	// Change the mesh
	_tessellation->RefineCells(PointsToRefine,directions,dr);
	// Fill the new hydro
	for(int i=N;i<N+n;++i)
	{
		_cells[i]=_cells[PointsToRefine[i-N]];
		CellsEvolve[i]=CellsEvolve[PointsToRefine[i-N]];
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
	vector<Vector2D> cor=_tessellation->GetMeshPoints();
	vector<int> order=HilbertOrder(cor,int(_cells.size()),innernum);
	ReArrangeVector(cor,order);
	if(cor.size()>order.size())
		cor.erase(cor.begin()+order.size(),cor.end());
	ReArrangeVector(_cells,order);
	if(tracer_flag_)
		ReArrangeVector(tracer_,order);
	_tessellation->Update(cor);
}
