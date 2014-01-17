#include <cmath>
#include "hydrodynamics_2d.hpp"
#include "../common/hydrodynamics.hpp"
#include "../../misc/universal_error.hpp"
#include <iostream>
#include <boost/scoped_ptr.hpp>

vector<Primitive> InitialiseCells
(SpatialDistribution const& density,
 SpatialDistribution const& pressure,
 SpatialDistribution const& xvelocity,
 SpatialDistribution const& yvelocity,
 EquationOfState const& eos,
 Tessellation const* tess,bool CMvalue)
{
  vector<Primitive> res(tess->GetPointNo());
  for(int i=0;i<(int)res.size();i++){
    Vector2D r;
    if(!CMvalue)
      r= tess->GetMeshPoint(i);
    else
      r = tess->GetCellCM(i);
    const Vector2D v(xvelocity.EvalAt(r),
		     yvelocity.EvalAt(r));
    res[i] = CalcPrimitive(density.EvalAt(r),
			   pressure.EvalAt(r),
			   v, eos);
  }
  return res;
}

vector<Conserved> CalcConservedIntensive
(vector<Primitive> const& cells)
{
  vector<Conserved> res(cells.size());
  for(int i=0;i<(int)cells.size();i++){
    res[i] = Primitive2Conserved(cells[i]);
  }
  return res;
}

vector<Conserved> CalcConservedExtensive
(vector<Conserved> const& cons_int,
 Tessellation const* tess)
{
  vector<Conserved> res(cons_int.size());
  for(int i=0;i<(int)cons_int.size();i++){
    res[i] = tess->GetVolume(i)*
      cons_int[i];
  }
  return res;
}

void CalcPointVelocities(Tessellation const* tessellation,
			 vector<Primitive> const& cells,
			 PointMotion *pointmotion,
			 vector<Vector2D>& pointvelocity,double time)
{
  /*for(int i = 0; i<tessellation->GetPointNo(); i++){    
    pointvelocity[i] = pointmotion->CalcVelocity(i, tessellation,cells,
    time);
    }*/
  pointvelocity=pointmotion->calcAllVelocities(tessellation,cells,time);
}


double TimeStepForCell(Primitive const& cell,double width,
		       vector<Vector2D> const& face_velocites)
{
  double Max=0;
  double temp;
  for(int i=0;i<(int)face_velocites.size();++i)
    {
      temp=abs(face_velocites[i]-cell.Velocity);
      if(temp>Max)
	Max=temp;
    }
  return width/(cell.SoundSpeed+Max);
}

double TimeStepForCellBoundary(Primitive const& cell,
			       vector<Primitive> const& cells,double width,
			       vector<Vector2D> const& face_velocites,Tessellation const* tess,
			       HydroBoundaryConditions const* hbc,int index,double time)
{
  double Max=0;
  double temp;
  vector<int> edge_index=tess->GetCellEdges(index);
  Edge edge;
  for(int i=0;i<(int)face_velocites.size();++i)
    {
      edge=tess->GetEdge(edge_index[i]);
      if(hbc->IsBoundary(edge,tess))
	temp=max(abs(face_velocites[i]-cell.Velocity),
		 abs(face_velocites[i]-hbc->GetBoundaryPrimitive(edge,tess,
								 cells,time).Velocity));
      else
	temp=abs(face_velocites[i]-cell.Velocity);
      if(temp>Max)
	Max=temp;
    }
  return width/(cell.SoundSpeed+Max);
}


double CalcTimeStep(Tessellation const* tessellation,
		    vector<Primitive> const& cells,
		    vector<Vector2D> const& facevelocity,
		    HydroBoundaryConditions const* hbc,double time,vector<CustomEvolution*>
		    const& evolve)
{
  vector<int> face_index=tessellation->GetCellEdges(0);
  vector<Vector2D> face_vel(face_index.size());
  bool first_time=true;
  double dt=0;
  double dt_temp=0;
  for(int i = 0;i<tessellation->GetPointNo(); i++)
    {
      if(!evolve.empty())
	if(evolve[i]!=0)
	  if(evolve[i]->TimeStepRelevant()==false)
	    continue;
      face_index=tessellation->GetCellEdges(i);
      face_vel.resize(face_index.size());
      for(int j=0;j<(int)face_index.size();++j)
	face_vel[j]=facevelocity[face_index[j]];
      if(!tessellation->NearBoundary(i))
	dt_temp=TimeStepForCell(cells[i],tessellation->GetWidth(i),face_vel);
      else
	dt_temp = TimeStepForCellBoundary(cells[i],cells,
					  tessellation->GetWidth(i),face_vel,tessellation,hbc,i,time);
      if(first_time)
	{
	  first_time=false;
	  dt=dt_temp;
	  continue;
	}
      if(dt_temp<dt)
	dt=dt_temp;
    }
  return dt;
}

void UpdateConservedExtensive
(Tessellation const* tessellation,
 vector<Conserved> const& fluxes,
 double dt, 
 vector<Conserved>& conserved_extensive,
 HydroBoundaryConditions const* boundaryconditions)
{
  for(int i=0;i<tessellation->GetPointNo();++i){
    vector<int> cell_edge_indices = tessellation->GetCellEdges(i);
    if(!boundaryconditions->IsGhostCell(i,tessellation)){
      for(int j=0;j<(int)cell_edge_indices.size();++j){
	int edge_index = cell_edge_indices[j];
	Edge edge = tessellation->GetEdge(edge_index);
	Conserved delta = dt*edge.GetLength()*fluxes[edge_index];
	if(i==edge.GetNeighbor(0))
	  conserved_extensive[i] -= delta;
	else if(i==edge.GetNeighbor(1))
	  conserved_extensive[i] += delta;
	else
	  {
	    UniversalError eo("Error in UpdateConservedExtensive: Cell and edge are not mutual neighbors");
	    eo.AddEntry("Edge number",edge_index);
	    eo.AddEntry("cell number",i);
	    throw eo;
	  }
      }
    }
  }
}

void MoveMeshPoints(vector<Vector2D> const& pointvelocity,
		    double dt, Tessellation* tessellation)
{
  vector<Vector2D> points;
  points.resize(tessellation->GetPointNo());
  for(int i=0;i<tessellation->GetPointNo(); i++){
    points[i] = tessellation->GetMeshPoint(i) + 
      dt*pointvelocity[i];
  }
  tessellation->Update(points);
}

void UpdateConservedIntensive(Tessellation const* tessellation,
			      vector<Conserved> const& conservedextensive,
			      vector<Conserved>& conservedintensive)
{
  for(int i = 0; i<tessellation->GetPointNo(); i++){
    conservedintensive[i] = conservedextensive[i]/
      tessellation->GetVolume(i);
  }
}

void UpdatePrimitives
(vector<Conserved> const& conservedintensive,
 EquationOfState const& eos,vector<Primitive>& cells,
 vector<CustomEvolution*> const& CellsEvolve,vector<Primitive> &old_cells,
 bool densityfloor,double densitymin,double pressuremin,Tessellation const*
 tess,double time,vector<vector<double> > const& tracers)
{
  Conserved temp;
  for(int i=0;i < (int)cells.size(); i++){
    try
      {
	if(CellsEvolve[i]==0)
	  {
	    temp=conservedintensive[i];
	    if(densityfloor)
	      {
		if(temp.Mass<densitymin)
		  {
		    temp.Mass=densitymin;
		    temp.Energy=temp.Mass*(eos.dp2e(temp.Mass,old_cells[i].Pressure*
						    temp.Mass/old_cells[i].Density)+
					   0.5*pow(abs(old_cells[i].Velocity),2.0));
		    temp.Momentum=temp.Mass*old_cells[i].Velocity;
		  }
		double p=eos.de2p(temp.Mass,(temp.Energy-
					     0.5*ScalarProd(temp.Momentum,temp.Momentum)/temp.Mass)/
				  temp.Mass);
		if(p<pressuremin)
		  temp.Energy=0.5*ScalarProd(temp.Momentum,temp.Momentum)/
		    temp.Mass+temp.Mass*eos.dp2e(temp.Mass,pressuremin);
		cells[i]=Conserved2Primitive(temp,eos);
	      }
	    else
	      cells[i] = Conserved2Primitive(temp,eos);
	  }
	else
	  cells[i]=CellsEvolve[i]->UpdatePrimitive(conservedintensive,
						   &eos,old_cells,i,tess,time,tracers);
      }
    catch(UniversalError& eo)
      {
	eo.AddEntry("UpdatePrimitives data starts here",0);
	eo.AddEntry("cell index",(double)i);
	throw;
      }
  }
}

namespace {
  Primitive RotatePrimitive(Vector2D const& normaldir,
			    Vector2D const& paraldir,
			    Primitive const& p)
  {
    Primitive res = p;
    res.Velocity.Set(Projection(p.Velocity,normaldir),
		     Projection(p.Velocity,paraldir));
    return res;
  }

  Conserved RotateFluxBack(Conserved const& c,
			   Vector2D const& normaldir,
			   Vector2D const& paraldir)
  {
    Conserved res = c;
    res.Momentum = c.Momentum.x*normaldir/abs(normaldir)+
      c.Momentum.y*paraldir/abs(paraldir);
    return res;
  }
}

Conserved FluxInBulk(Vector2D const& normaldir,
		     Vector2D const& paraldir,
		     Primitive const& left,
		     Primitive const& right,
		     Vector2D const& edge_velocity,
		     RiemannSolver const& rs)
{
  Primitive rotated_left = RotatePrimitive(normaldir, paraldir, left);
  Primitive rotated_right = RotatePrimitive(normaldir, paraldir, right);
  double normal_speed = Projection(edge_velocity,normaldir);
  Conserved res = rs.Solve(rotated_left, rotated_right, normal_speed);
  return RotateFluxBack(res, normaldir, paraldir);
}

namespace {
  vector<bool> is_normal_evolution(vector<CustomEvolution*> cev)
  {
    vector<bool> res(cev.size(),true);
    /*for(int i=0;i<(int)cev.size();++i)
      res[i] = cev[i]==0;*/
    return res;
  }
}

namespace {
  bool primitive_is_nan(Primitive const& p)
  {
    return is_nan(p.Density)||
      is_nan(p.Pressure)||
      is_nan(p.Velocity.x)||
      is_nan(p.Velocity.y)||
      is_nan(p.Energy)||
      is_nan(p.SoundSpeed);
  }
}

void CalcFluxes
(Tessellation const* tessellation,
 vector<Primitive> const& cells,
 double dt,
 double time,
 SpatialReconstruction* interpolation,
 vector<Vector2D> const& facevelocity,
 HydroBoundaryConditions const* boundaryconditions,
 RiemannSolver const& rs,
 vector<Conserved>& fluxes,
 vector<CustomEvolution*> const& CellsEvolve,
 vector<vector<double> > const& tracers)
{
  fluxes.resize(tessellation->GetTotalSidesNumber());

  try
    {
      interpolation->Prepare(tessellation,cells,tracers,dt,time);
    }
  catch(UniversalError& eo)
    {
      eo.AddEntry("Error in CalcFlux Interpolation prepare",0);
      throw;
    }
  for(int i = 0; i < tessellation->GetTotalSidesNumber(); i++)
    {
      try
	{
	  const Edge edge = tessellation->GetEdge(i);
	  const int n1=edge.GetNeighbor(1);
	  const int n0=edge.GetNeighbor(0);
	  if(!boundaryconditions->IsBoundary(edge,tessellation))
	    {
	      if(CellsEvolve[n0]==0&&CellsEvolve[n1]==0)
		{
		  // In the bulk of the grid

		  const Vector2D normaldir = 
		    tessellation->GetMeshPoint(edge.GetNeighbor(1))-
		    tessellation->GetMeshPoint(edge.GetNeighbor(0));

		  const Vector2D paraldir = 
		    edge.GetVertex(1) - edge.GetVertex(0);
			      
		  const Primitive left = interpolation->Interpolate
		    (tessellation, cells, dt, edge, 0,InBulk,facevelocity[i]);
		  const Primitive right = interpolation->Interpolate
		    (tessellation, cells, dt, edge, 1,InBulk,facevelocity[i]);
		  if(primitive_is_nan(right)||primitive_is_nan(left)){
		    UniversalError eo("Nan occurred during flux calculation right after interpolation");
		    throw;
		  }
		  fluxes[i] = FluxInBulk(normaldir, paraldir,
					 left, right,
					 facevelocity[i], rs);
		}
	      else
		{
		  if(CellsEvolve[n0]!=0&&CellsEvolve[n1]!=0&&
		     CellsEvolve[n0]!=CellsEvolve[n1]&&
		     (!CellsEvolve[n0]->flux_indifferent()||
		      !CellsEvolve[n1]->flux_indifferent()))
		    {
		      UniversalError eo("Bad inner boundary, edge doesn't have same cellevolve conditions on both sides");
		      eo.AddEntry("First cell",n0);
		      eo.AddEntry("Second cell",n1);
		      eo.AddEntry("Edge number",i);
		      throw eo;
		    }
		  if(CellsEvolve[n0]!=0)
		    fluxes[i]=CellsEvolve[n0]->CalcFlux
		      (tessellation,cells,
		       dt,interpolation,edge,facevelocity[i],
		       rs,n0,boundaryconditions,time,tracers);
		  else
		    fluxes[i]=CellsEvolve[n1]->CalcFlux
		      (tessellation,cells,
		       dt,interpolation,edge,facevelocity[i],
		       rs,n1,boundaryconditions,time,tracers);
		}
	    }
	  else
	    {
	      // Boundaries
	      fluxes[i] = boundaryconditions->CalcFlux
		(tessellation, cells, facevelocity[i], edge,
		 interpolation,dt,time);
	    }
	}
      catch(UniversalError& eo)
	{
	  eo.AddEntry("Error in CalcFlux",0);
	  eo.AddEntry("edge index",(double)i);
	  eo.AddEntry("edge x1 location = ",tessellation->GetEdge(i).get_x(0));
	  eo.AddEntry("edge y1 location = ",tessellation->GetEdge(i).get_y(0));
	  eo.AddEntry("edge x2 location = ",tessellation->GetEdge(i).get_x(1));
	  eo.AddEntry("edge y2 location = ",tessellation->GetEdge(i).get_y(1));
	  throw;
	}

    }
}

void ExternalForceContribution
(Tessellation const* tess,
 vector<Primitive> const& cells,SourceTerm *force,
 double t,
 double dt,
 vector<Conserved>& conserved_extensive,
 HydroBoundaryConditions const* hbc,vector<Conserved> const& fluxes,
 vector<Vector2D> const& point_velocity,vector<double> &g,bool coldflows_flag,
 vector<vector<double> > &tracers_extensive)
{
  vector<double> dtracer;
  if(!tracers_extensive.empty())
    dtracer.assign(tracers_extensive[0].size(),0);
  if(coldflows_flag)
    g.resize(point_velocity.size());
  for(int i=0;i<tess->GetPointNo();++i)
    {
      if(!hbc->IsGhostCell(i,tess))
	{
	  Conserved cons(force->Calculate
			 (tess,cells,i,fluxes,point_velocity,hbc,tracers_extensive,
			  dtracer,t,dt));
	  conserved_extensive[i]+=dt*cons;
	  if(!tracers_extensive.empty())
	    {
	      tracers_extensive[i]=tracers_extensive[i]+dt*dtracer;
	    }
	  if(coldflows_flag)
	    g[i]=abs(cons.Momentum)/(cells[i].Density*tess->GetVolume(i));
	}
    }
}


HydroSnapshot::HydroSnapshot
(vector<Vector2D> const& mesh_points_i,
 vector<Primitive> const& cells_i):
  mesh_points(mesh_points_i),
  cells(cells_i) {}

HydroSnapshot::HydroSnapshot(void):
  mesh_points(vector<Vector2D>()),
  cells(vector<Primitive>()) {}

vector<Vector2D> calc_point_velocities
(Tessellation const* tess,
 vector<Primitive> const& cells,
 PointMotion *point_motion,
 double time)
{
  return point_motion->calcAllVelocities(tess,cells,time);
}

namespace {
  vector<Conserved> calc_fluxes
  (Tessellation const* tess,
   vector<Primitive> const& cells,
   double dt,
   double time,
   SpatialReconstruction* interpolation,
   vector<Vector2D> const& edge_velocities,
   HydroBoundaryConditions const* hbc,
   RiemannSolver const& rs,
   vector<CustomEvolution*> const& CellsEvolve,
   vector<vector<double> > const& tracers)
  {
    vector<Conserved> res;
    CalcFluxes(tess, cells, dt,
	       time,
	       interpolation,
	       edge_velocities,
	       hbc, rs, res,CellsEvolve,tracers);
    return res;
  }
}
vector<Vector2D> get_all_mesh_points
(Tessellation const* tess)
{
  vector<Vector2D> res(tess->GetPointNo());
  for(int i=0;i<(int)tess->GetPointNo();++i)
    res[i] = tess->GetMeshPoint(i);
  return res;
}

double TimeAdvance2mid
(Tessellation* tess,vector<Primitive> &cells,
 PointMotion *point_motion,HydroBoundaryConditions const* hbc,
 SpatialReconstruction* interpolation,RiemannSolver const& rs,
 EquationOfState const& eos,SourceTerm* force,double time,double cfl,
 double endtime,vector<CustomEvolution*> const& CellsEvolve,
 vector<vector<double> > &tracers,
 double dt_external,bool traceflag,bool coldflows_flag,
 double as,double bs,bool densityfloor,double densitymin,
 double pressuremin,bool EntropyCalc)
{
  boost::scoped_ptr<Tessellation> sp_tess(tess->clone());
  Tessellation *tess_temp = sp_tess.get();

  //do half time step
  vector<Vector2D> point_velocities = calc_point_velocities
    (tess_temp, cells, point_motion, time);

  vector<Vector2D> edge_velocities =tess_temp->calc_edge_velocities
    (hbc,
     point_velocities,time);

  double dt = cfl*CalcTimeStep
    (tess_temp,cells, edge_velocities,hbc,time,CellsEvolve);

  if(dt_external>0)
    dt=min(dt,dt_external*cfl);

  if(endtime>0)
    if(time+dt>endtime)
      dt=endtime-time;

  // Entropy and tracers evolution
  vector<double> g,Ek,Ef;
  if(coldflows_flag)
    {
      int n=int(cells.size());
      for(int i=0;i<n;++i)
	tracers[i][0]=eos.dp2s(cells[i].Density,cells[i].Pressure);
    }

  vector<Conserved> fluxes = calc_fluxes
    (tess_temp, cells, 0.5*dt, time, interpolation,
     edge_velocities, hbc, rs,CellsEvolve,tracers);

  vector<vector<double> > old_trace=tracers;
  vector<vector<double> > tracer_extensive;
  if(traceflag)
    {
      vector<vector<double> > trace_change;
      trace_change = CalcTraceChange(tracers,cells,tess_temp,fluxes,0.5*dt,hbc,
				     interpolation,time,CellsEvolve,edge_velocities);
      MakeTracerExtensive(tracers,tess_temp,cells,tracer_extensive);
      UpdateTracerExtensive(tracer_extensive,trace_change,CellsEvolve,cells,
			    tess_temp,time);
    }

  vector<Conserved> intensive = CalcConservedIntensive(cells);

  vector<Conserved> extensive = CalcConservedExtensive
    (intensive, tess_temp);

  UpdateConservedExtensive(tess_temp, fluxes, 0.5*dt,
			   extensive, hbc);

  ExternalForceContribution
    (tess_temp, cells,force, time, 0.5*dt,
     extensive, hbc,fluxes,point_velocities,g,coldflows_flag,tracers);

  MoveMeshPoints(point_velocities, 0.5*dt, tess_temp);

  UpdateConservedIntensive(tess_temp, extensive, intensive);


  if(coldflows_flag)
    {
      Ek=GetMaxKineticEnergy(tess_temp,cells,CellsEvolve);
      Ef=GetForceEnergy(tess_temp,g);
      FixPressure(intensive,tracer_extensive,eos,Ek,Ef,as,bs,CellsEvolve,
		  tess_temp,extensive,cells,hbc,time);
    }

  vector<Primitive> new_cells(cells.size());
  UpdatePrimitives(intensive, eos, new_cells,CellsEvolve,cells,densityfloor,
		   densitymin,pressuremin,tess_temp,time+0.5*dt,tracers);
  if(traceflag)
    {
      MakeTracerIntensive(tracers,tracer_extensive,tess_temp,new_cells);
    }

  // End half step

  point_velocities = calc_point_velocities(tess_temp,new_cells,
					   point_motion, time+0.5*dt);
  edge_velocities = tess_temp->calc_edge_velocities(hbc,point_velocities,time);

  fluxes = calc_fluxes
    (tess_temp, new_cells, dt, time, interpolation,
     edge_velocities, hbc, rs,CellsEvolve,tracers);

  if(coldflows_flag&&EntropyCalc)
    {
      int n=int(cells.size());
      for(int i=0;i<n;++i)
	tracers[i][0]=eos.dp2s(new_cells[i].Density,new_cells[i].Pressure);
    }
  if(traceflag)
    {
      vector<vector<double> > trace_change;
      trace_change = CalcTraceChange(tracers,new_cells,tess_temp,fluxes,dt,hbc,
				     interpolation,time,CellsEvolve,edge_velocities);
      tracers=old_trace;
      MakeTracerExtensive(tracers,tess,cells,tracer_extensive);
      UpdateTracerExtensive(tracer_extensive,trace_change,CellsEvolve,new_cells,
			    tess_temp,time);
    }

  intensive = CalcConservedIntensive(cells);

  extensive = CalcConservedExtensive
    (intensive, tess);

  UpdateConservedExtensive(tess_temp, fluxes, dt,
			   extensive, hbc);

  ExternalForceContribution
    (tess_temp, new_cells,force, time+0.5*dt,dt,
     extensive, hbc,fluxes,point_velocities,g,coldflows_flag,tracers);

  //Move mesh points 
  MoveMeshPoints(point_velocities,dt,tess);

  UpdateConservedIntensive(tess, extensive, intensive);


  if(coldflows_flag)
    {
      Ek=GetMaxKineticEnergy(tess,cells,CellsEvolve);
      Ef=GetForceEnergy(tess,g);
      FixPressure(intensive,tracer_extensive,eos,Ek,Ef,as,bs,CellsEvolve,tess,
		  extensive,cells,hbc,time);
    }

  UpdatePrimitives
    (intensive, eos, cells,CellsEvolve,cells,densityfloor,
     densitymin,pressuremin,tess,time+dt,tracers);

  if(traceflag)
    {
      MakeTracerIntensive(tracers,tracer_extensive,tess,cells);
    }


  //	delete tess_temp;
  //tess_temp.reset();
  return dt;
}

vector<Primitive> make_eos_consistent
(vector<Primitive> const& vp,
 EquationOfState const& eos)
{
  vector<Primitive> res = vp;
  for(int i=0;i<(int)vp.size();++i)
    res[i] = make_eos_consistent(vp[i],eos);
  return res;
}

vector<vector<double> > CalcTraceChange
(vector<vector<double> > const& old_trace,vector<Primitive> const& cells,
 Tessellation const* tess,vector<Conserved> const& fluxes,double dt,
 HydroBoundaryConditions const* hbc, SpatialReconstruction const* interp,
 double time,vector<CustomEvolution*> const& CellsEvolve,
 vector<Vector2D> const& edge_velocities)
{
  vector<vector<double> > res;
  int n=int(old_trace.size());
  int dim;
  if(n>0)
    dim=int(old_trace[0].size());
  else
    return res;
  res.resize(n);
  double dm;
  vector<int> cell_edge_indices;
  Edge edge;
  int edge_index;
  vector<double> temp(dim);
  // calculate change
  for(int i=0;i<n;++i)
    {
      res[i].resize(dim,0);
      cell_edge_indices=tess->GetCellEdges(i);
      int Nedges=int(cell_edge_indices.size());
      int n1,n0;
      for(int k=0;k<Nedges;++k)
	{
	  edge_index = cell_edge_indices[k];
	  edge = tess->GetEdge(edge_index);
	  dm=fluxes[edge_index].Mass;
	  n1=edge.GetNeighbor(1);
	  n0=edge.GetNeighbor(0);
	  if(!hbc->IsBoundary(edge,tess))
	    {
	      if(CellsEvolve[n0]!=0&&CellsEvolve[n1]!=0&&
		 CellsEvolve[n0]!=CellsEvolve[n1]&&
		 (!CellsEvolve[n0]->flux_indifferent()||
		  !CellsEvolve[n1]->flux_indifferent()))
		{
		  UniversalError eo("Bad inner boundary, edge doesn't have same cellevolve conditions on both sides");
		  eo.AddEntry("First cell",n0);
		  eo.AddEntry("Second cell",n1);
		  throw eo;
		}
	      if(CellsEvolve[n0]!=0)
		{
		  temp=CellsEvolve[n0]->CalcTracerFlux(tess,cells,old_trace,dm,edge,n0,dt,
						       time,interp,edge_velocities[edge_index]);
		  transform(res[i].begin(),res[i].end(),temp.begin(),res[i].begin(),
			    plus<double>());
		  continue;
		}
	      if(CellsEvolve[n1]!=0)
		{
		  temp=CellsEvolve[n1]->CalcTracerFlux(tess,cells,old_trace,dm,edge,n1,dt,
						       time,interp,edge_velocities[edge_index]);
		  transform(res[i].begin(),res[i].end(),temp.begin(),res[i].begin(),
			    plus<double>());
		  continue;
		}
	    }
	  if(dm>0)
	    {
	      if(i==n1)
		{
		  if(hbc->IsBoundary(edge,tess))
		    {
		      temp=hbc->CalcTracerFlux(tess,cells,old_trace,dm,edge,i,dt,time,interp,
					       edge_velocities[edge_index]);
		      transform(res[i].begin(),res[i].end(),temp.begin(),res[i].begin(),
				plus<double>());
		    }
		  else
		    {
		      temp=interp->interpolateTracers(tess,cells,old_trace,dt,edge,0,InBulk,
						      edge_velocities[edge_index]);
		      transform(temp.begin(),temp.end(),temp.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		      transform(res[i].begin(),res[i].end(),temp.begin(),res[i].begin(),
				plus<double>());
		    }
		}
	      else
		{
		  if(hbc->IsBoundary(edge,tess))
		    {
		      temp=hbc->CalcTracerFlux(tess,cells,old_trace,dm,edge,i,dt,time,interp,
					       edge_velocities[edge_index]);
		      transform(res[i].begin(),res[i].end(),temp.begin(),res[i].begin(),
				minus<double>());
		    }
		  else
		    {
		      temp=interp->interpolateTracers(tess,cells,old_trace,dt,edge,
						      0,InBulk,edge_velocities[edge_index]);
		      transform(temp.begin(),temp.end(),temp.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		      transform(res[i].begin(),res[i].end(),temp.begin(),res[i].begin(),
				minus<double>());
		    }
		}
	    }
	  else
	    {
	      if(i==n1)
		{
		  if(hbc->IsBoundary(edge,tess))
		    {
		      temp=hbc->CalcTracerFlux(tess,cells,old_trace,dm,edge,i,dt,
					       time,interp,edge_velocities[edge_index]);
		      transform(res[i].begin(),res[i].end(),temp.begin(),res[i].begin(),
				plus<double>());
		    }
		  else
		    {
		      temp=interp->interpolateTracers(tess,cells,old_trace,dt,
						      edge,1,InBulk,edge_velocities[edge_index]);
		      transform(temp.begin(),temp.end(),temp.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		      transform(res[i].begin(),res[i].end(),temp.begin(),res[i].begin(),
				plus<double>());
		    }
		}
	      else
		{
		  if(hbc->IsBoundary(edge,tess))
		    {
		      temp=hbc->CalcTracerFlux(tess,cells,old_trace,dm,edge,i,dt,
					       time,interp,edge_velocities[edge_index]);
		      transform(res[i].begin(),res[i].end(),temp.begin(),res[i].begin(),
				minus<double>());
		    }
		  else
		    {
		      temp=interp->interpolateTracers(tess,cells,old_trace,dt,
						      edge,1,InBulk,edge_velocities[edge_index]);
		      transform(temp.begin(),temp.end(),temp.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		      transform(res[i].begin(),res[i].end(),temp.begin(),res[i].begin(),
				minus<double>());
		    }
		}
	    }
	}
    }
  return res;
}

vector<double> GetMaxKineticEnergy(Tessellation const* tess,vector<Primitive> const&
				   cells,vector<CustomEvolution*> const& customevolve)
{
  int n=int(cells.size());
  vector<double> res;
  res.resize(n);
  double e;
  for(int j=0;j<n;++j)
    {
      e=0;
      if(!NearBoundary(j,tess,customevolve))
	{	
	  vector<int> neigh=tess->GetNeighbors(j);
	  e=pow(abs(cells[j].Velocity-cells[neigh[0]].Velocity),2);
	  for(int i=1;i<(int)neigh.size();++i)
	    {// This could be made much faster by writing the expression implicitly
	      e=max(e,pow(abs(cells[j].Velocity-cells[neigh[i]].Velocity),2));
	    }
	}
      res[j]=0.5*e;
    }
  return res;
}

vector<double> GetForceEnergy(Tessellation const* tess,
			      vector<double> const& g)
{
  vector<double> res;
  int n=int(g.size());
  res.resize(n);
  for(int i=0;i<n;++i)
    res[i]=g[i]*tess->GetWidth(i);
  return res;
}

void FixPressure(vector<Conserved> &intensive,vector<vector<double> > const& entropy,
		 EquationOfState const& eos,vector<double> const& Ek,
		 vector<double> const& Ef,double as,double bs,vector<CustomEvolution*>
		 const& customevolve,Tessellation const* tess,vector<Conserved> &extensive,
		 vector<Primitive> const& cells,HydroBoundaryConditions const* hbc,
		 double time)
{
  int n=int(intensive.size());
  double Et,Ek2;
  double temp;
  for(int i=0;i<n;++i){
    if(customevolve[i]==0){
      //Make intensive
      temp=entropy[i][0]/(tess->GetVolume(i)*intensive[i].Mass);
      Ek2=0.5*pow(abs(intensive[i].Momentum)/intensive[i].Mass,2);
      Et=intensive[i].Energy/intensive[i].Mass-Ek2;
      if((Et<as*Ek[i])||(Et<bs*Ef[i])){
	if(~IsShockedCell(tess,i,cells,hbc,time)||Et<0){
	  Et=eos.dp2e(intensive[i].Mass,
		      eos.sd2p(temp,intensive[i].Mass));
	  intensive[i].Energy=intensive[i].Mass*(Et+Ek2);
	  extensive[i].Energy=tess->GetVolume(i)*intensive[i].Energy;
	}
      }
    }
  }
}

bool NearBoundary(int index,Tessellation const* tess,
		  vector<CustomEvolution*> const& /*customevolve*/)
{
  vector<int> neigh=tess->GetNeighbors(index);
  int n=int(neigh.size());
  for(int i=0;i<n;++i)
    {
      if(neigh[i]<0)
	return true;
      /*if(customevolve[neigh[i]]!=0)
	return true;*/
    }
  return false;
}

void MakeTracerExtensive(vector<vector<double> > const &tracer,
			 Tessellation const* tess,vector<Primitive> const& cells,vector<vector<double> >
			 &result)
{
  int n=int(tracer.size());
  if(n==0)
    return;
  int dim=int(tracer[0].size());
  double mass;
  result.resize(n);
  for(int i=0;i<n;++i)
    {
      mass=tess->GetVolume(i)*cells[i].Density;
      result[i].resize(dim);
      for(int j=0;j<dim;++j)
	result[i][j]=tracer[i][j]*mass;
    }
  return;
}

void MakeTracerIntensive(vector<vector<double> > &tracer,vector<vector<double> >
			 const& tracer_extensive,Tessellation const* tess,vector<Primitive> const& cells)
{
  int n=int(tracer.size());
  if(n==0)
    return;
  int dim=int(tracer[0].size());
  double mass;
  for(int i=0;i<n;++i)
    {
      mass=tess->GetVolume(i)*cells[i].Density;
      for(int j=0;j<dim;++j)
	tracer[i][j]=tracer_extensive[i][j]/mass;
    }
  return;
}

void UpdateTracerExtensive(vector<vector<double> > &tracerextensive,
			   vector<vector<double> > const& tracerchange,vector<CustomEvolution*> const&
			   CellsEvolve,vector<Primitive> const& cells,Tessellation const* tess,
			   double time)
{
  int n=int(tracerextensive.size());
  if(n==0)
    return;
  int dim=int(tracerextensive[0].size());
  for(int i=0;i<n;++i)
    if(CellsEvolve[i]!=0)
      {
	tracerextensive[i]=CellsEvolve[i]->UpdateTracer(
							i,tracerextensive,cells,tess,time);
      }
    else
      for(int j=0;j<dim;++j)
	tracerextensive[i][j]+=tracerchange[i][j];
  return;
}

void TracerResetCalc(double alpha,SpatialDistribution const& originalD,
		     SpatialDistribution const& originalP,SpatialDistribution const& originalVx,
		     SpatialDistribution const& originalVy, vector<Primitive> &cells,
		     Tessellation const* tess,vector<vector<double> > &tracer,
		     int tracerindex,EquationOfState const& eos,int InnerNum)
{
  int n=int(cells.size());
  if(n<1)
    return;
  Vector2D velocity;
  if(tracer.empty())
    {
      for(int i=0;i<InnerNum;++i)
	{
	  velocity.Set(originalVx.EvalAt(tess->GetCellCM(i)),
		       originalVy.EvalAt(tess->GetCellCM(i)));
	  cells[i]=CalcPrimitive(originalD.EvalAt(tess->GetCellCM(i)),
				 originalP.EvalAt(tess->GetCellCM(i)),velocity,eos);
	}
      return;	
    }
  if(tracerindex>=(int)tracer[0].size()||tracerindex<0)
    throw UniversalError("Error in tracerReset, wrong dimension for tracer");
  for(int i=0;i<n;++i)
    {
      if((tracer[i][tracerindex]<alpha)||(i<InnerNum))//*cells[i].Density)
	{
	  velocity.Set(originalVx.EvalAt(tess->GetCellCM(i)),
		       originalVy.EvalAt(tess->GetCellCM(i)));
	  cells[i]=CalcPrimitive(originalD.EvalAt(tess->GetCellCM(i)),
				 originalP.EvalAt(tess->GetCellCM(i)),velocity,eos);
	  if(tracer[i][tracerindex]<0)
	    tracer[i][tracerindex]=0;
	  if(i<InnerNum)
	    tracer[i][tracerindex]=0;
	}
    }
  return;
}

void GetPointToRemove(Tessellation const* tess,Vector2D const& point,
		      double R,vector<int> & PointToRemove,int Inner)
{
  int n=tess->GetPointNo();
  PointToRemove.clear();
  bool test;
  for(int i=Inner;i<n;++i)
    {
      // Check if point is completly engulfed
      test=true;
      vector<int> neigh=tess->GetNeighbors(i);
      for(int j=0;j<(int)neigh.size();++j)
	if(neigh[j]>=Inner)
	  test=false;
      // Is point inside a radius?
      if(abs(point-tess->GetMeshPoint(i))<R||test)
	PointToRemove.push_back(i);
    }
  return;
}

namespace {
  Vector2D GetReflectedPoint(Tessellation const* tess,int point,
			     Edge const& edge)
  {
    Vector2D MeshPoint=tess->GetMeshPoint(point);
    Vector2D par=edge.GetVertex(1)-edge.GetVertex(0);
    if(abs(par.x)>abs(par.y)){
      // We are in the x direction
      if(MeshPoint.y>edge.get_y(0))
	MeshPoint.y-=MeshPoint.y-edge.get_y(0);
      else
	MeshPoint.y+=edge.get_y(0)-MeshPoint.y;
    }
    else{
      // We are in the y direction
      if(MeshPoint.x>edge.get_x(0))
	MeshPoint.x-=MeshPoint.x-edge.get_x(0);
      else
	MeshPoint.x+=edge.get_x(0)-MeshPoint.x;
    }
    return MeshPoint;
  }
}

bool IsShockedCell(Tessellation const* tess,int index,
		   vector<Primitive> const& cells,HydroBoundaryConditions const* hbc,
		   double time)
{
  double DivV=0;
  vector<int> edges_loc=tess->GetCellEdges(index);
  vector<Edge> edges(edges_loc.size());
  for(int i=0;i<(int)edges.size();++i)
    edges[i]=tess->GetEdge(edges_loc[i]);

  // Calculate gradiant by gauss theorem
  double vx_i = cells[index].Velocity.x;
  double vy_i = cells[index].Velocity.y;
  int other;
  Vector2D center=tess->GetMeshPoint(index);
  Vector2D c_ij,r_ij;
  for(int i=0;i<(int)edges.size();i++)
    {
      if(hbc->IsBoundary(edges[i],tess))
	{
	  Primitive other2=hbc->GetBoundaryPrimitive(edges[i],tess,cells,
						     time);
	  if((edges[i].GetNeighbor(0)==-1)||(edges[i].GetNeighbor(1)==-1))
	    {
	      c_ij = CalcCentroid(edges[i])-0.5*(GetReflectedPoint
						 (tess,index,edges[i])+center);
	      r_ij = center-GetReflectedPoint(tess,index,edges[i]);
	    }
	  else
	    {
	      int other_point;
	      if(edges[i].GetNeighbor(0)==index)
		other_point=edges[i].GetNeighbor(1);
	      else
		other_point=edges[i].GetNeighbor(0);
	      c_ij = CalcCentroid(edges[i])-0.5*(
						 tess->GetMeshPoint(other_point)+center);
	      r_ij = center-tess->GetMeshPoint(other_point);
	    }
	  double rij_1=1/abs(r_ij);
	  double vx_j = other2.Velocity.x;
	  double vy_j = other2.Velocity.y;
	  DivV+=edges[i].GetLength()*((vx_j-vx_i)*c_ij.x-0.5*
				      (vx_i+vx_j)*r_ij.x+(vy_j-vy_i)*c_ij.y-0.5*
				      (vy_i+vy_j)*r_ij.y)*rij_1;
	}
      else
	{
	  if(edges[i].GetNeighbor(0)==index)
	    other=edges[i].GetNeighbor(1);
	  else
	    other=edges[i].GetNeighbor(0);	
	  c_ij = CalcCentroid(edges[i])-
	    0.5*(tess->GetMeshPoint(other)+center);
	  r_ij = center-tess->GetMeshPoint(other);
	  double rij_1=1/abs(r_ij);
	  double vx_j = cells[other].Velocity.x;
	  double vy_j = cells[other].Velocity.y;
	  DivV+=edges[i].GetLength()*((vx_j-vx_i)*c_ij.x-0.5*
				      (vx_i+vx_j)*r_ij.x+(vy_j-vy_i)*c_ij.y-0.5*
				      (vy_i+vy_j)*r_ij.y)*rij_1;
	}
    }
  if(DivV<0)
    return true;
  else
    return false;
}
