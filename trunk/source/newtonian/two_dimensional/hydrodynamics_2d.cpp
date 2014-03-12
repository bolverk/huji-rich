#include <cmath>
#include "hydrodynamics_2d.hpp"
#include "../common/hydrodynamics.hpp"
#include "../../misc/universal_error.hpp"
#include "../../misc/utils.hpp"
#include <boost/scoped_ptr.hpp>

namespace {
  Vector2D calc_representing_point(Tessellation const& tess,
				   int index,
				   bool cm_flag)
  {
    if(cm_flag)
      return tess.GetCellCM(index);
    else
      return tess.GetMeshPoint(index);
  }

  Primitive initialize_single_cell(Tessellation const& tess,
				   int index,
				   bool cm_flag,
				   SpatialDistribution const& density,
				   SpatialDistribution const& pressure,
				   SpatialDistribution const& xvelocity,
				   SpatialDistribution const& yvelocity,
				   EquationOfState const& eos)
  {
    const Vector2D r = calc_representing_point(tess,index,cm_flag);
    return CalcPrimitive(density.EvalAt(r),
			 pressure.EvalAt(r),
			 Vector2D(xvelocity.EvalAt(r),
				  yvelocity.EvalAt(r)),
			 eos);
  }

  class CellInitializer: public VectorTermInitializer<Primitive>
  {
  public:

    CellInitializer(Tessellation const& tess,
		    bool cm_flag,
		    SpatialDistribution const& density,
		    SpatialDistribution const& pressure,
		    SpatialDistribution const& xvelocity,
		    SpatialDistribution const& yvelocity,
		    EquationOfState const& eos):
      tess_(tess), cm_flag_(cm_flag),
      density_(density),
      pressure_(pressure),
      xvelocity_(xvelocity),
      yvelocity_(yvelocity),
      eos_(eos) {}

    Primitive operator()(int n) const
    {
      return initialize_single_cell(tess_,
				    n,
				    cm_flag_,
				    density_,
				    pressure_,
				    xvelocity_,
				    yvelocity_,
				    eos_);
    }

    int getLength(void) const
    {
      return tess_.GetPointNo();
    }

    ~CellInitializer(void) {}

  private:
    Tessellation const& tess_;
    const bool cm_flag_;
    SpatialDistribution const& density_;
    SpatialDistribution const& pressure_;
    SpatialDistribution const& xvelocity_;
    SpatialDistribution const& yvelocity_;
    EquationOfState const& eos_;
  };
}

vector<Primitive> InitialiseCells
(SpatialDistribution const& density,
 SpatialDistribution const& pressure,
 SpatialDistribution const& xvelocity,
 SpatialDistribution const& yvelocity,
 EquationOfState const& eos,
 Tessellation const& tess,
 bool cm_value)
{
  return initialize_vector<Primitive>
    (CellInitializer(tess,
		     cm_value,
		     density,
		     pressure,
		     xvelocity,
		     yvelocity,
		     eos));		     
}

namespace {
  class IntensiveInitializer: public VectorTermInitializer<Conserved>
  {
  public:

    IntensiveInitializer(vector<Primitive> const& cells):
      cells_(cells) {}

    Conserved operator()(int n) const
    {
      return Primitive2Conserved(cells_[n]);
    }

    int getLength(void) const
    {
      return (int)cells_.size();
    }

  private:
    vector<Primitive> const& cells_;
  };
}

vector<Conserved> CalcConservedIntensive
(vector<Primitive> const& cells)
{
  return initialize_vector<Conserved>
    (IntensiveInitializer(cells));
}

namespace {
  class ExtensiveInitializer: public VectorTermInitializer<Conserved>
  {
  public:

    ExtensiveInitializer(vector<Conserved> const& intensive,
			 Tessellation const& tess):
      intensive_(intensive), tess_(tess) {}

    Conserved operator()(int n) const
    {
      return tess_.GetVolume(n)*intensive_[n];
    }

    int getLength(void) const
    {
      return tess_.GetPointNo();
    }

  private:
    vector<Conserved> const& intensive_;
    Tessellation const& tess_;
  };
}

vector<Conserved> CalcConservedExtensive
(vector<Conserved> const& cons_int,
 Tessellation const& tess)
{
  return initialize_vector<Conserved>
    (ExtensiveInitializer(cons_int, tess));
}

void CalcPointVelocities(Tessellation const& tessellation,
			 vector<Primitive> const& cells,
			 PointMotion& pointmotion,
			 vector<Vector2D>& pointvelocity,double time)
{
  pointvelocity=pointmotion.calcAllVelocities(tessellation,cells,time);
}


double TimeStepForCell(Primitive const& cell,double width,
		       vector<Vector2D> const& face_velocites)
{
  double max_fv=0;
  for(int i=0;i<(int)face_velocites.size();++i)
    max_fv = max(max_fv,abs(face_velocites[i]-cell.Velocity));
  return width/(cell.SoundSpeed+max_fv);
}

namespace {
  double calc_boundary_face_velocity(HydroBoundaryConditions const& hbc,
				     Edge const& edge,
				     Tessellation const& tess,
				     vector<Vector2D> const& face_velocities,
				     Primitive const& cell,
				     vector<Primitive> const& cells,
				     double time,
				     int i)
  {
    if(hbc.IsBoundary(edge,tess))
      return max(abs(face_velocities[i]-cell.Velocity),
		 abs(face_velocities[i]-
		     hbc.GetBoundaryPrimitive(edge,
					      tess,
					      cells,
					      time).Velocity));
    else
      return abs(face_velocities[i]-cell.Velocity);
  }
}

double TimeStepForCellBoundary
(Primitive const& cell,
 vector<Primitive> const& cells,
 double width,
 vector<Vector2D> const& face_velocities,
 Tessellation const& tess,
 HydroBoundaryConditions const& hbc,
 int index,double time)
{
  double max_fv=0;
  const vector<int> edge_index=tess.GetCellEdges(index);
  for(int i=0;i<(int)face_velocities.size();++i)
    {
      const Edge edge=tess.GetEdge(edge_index[i]);
      const double temp = calc_boundary_face_velocity
	(hbc, edge, tess, face_velocities, cell, 
	 cells, time, i);
      max_fv = max(max_fv,temp);
    }
  return width/(cell.SoundSpeed+max_fv);
}

namespace {
  bool irrelevant_for_time_step(vector<CustomEvolution*> const& cevolve,
				int i)
  {
    if(!cevolve.empty())
      if(cevolve[i]!=0)
	if(!cevolve[i]->TimeStepRelevant())
	  return true;
    return false;
  }
}

namespace {
  double calc_dt_temp(Tessellation const& tess,
		      vector<Primitive> const& cells,
		      vector<Vector2D> const& face_vel,
		      HydroBoundaryConditions const& hbc,
		      double time,
		      int i)
  {
    if(!tess.NearBoundary(i))
      return TimeStepForCell(cells[i],tess.GetWidth(i),face_vel);
    else
      return TimeStepForCellBoundary(cells[i],
				     cells,
				     tess.GetWidth(i),
				     face_vel,
				     tess,
				     hbc,
				     i,
				     time);
  }

  class FaceVelocityInitializer: public VectorTermInitializer<Vector2D>
  {
  public:

    FaceVelocityInitializer(vector<int> const& face_index,
			    vector<Vector2D> const& face_velocity):
      face_index_(face_index),
      face_velocity_(face_velocity) {}

    Vector2D operator()(int n) const
    {
      return face_velocity_[face_index_[n]];
    }

    int getLength(void) const
    {
      return (int)face_index_.size();
    }

  private:
    vector<int> const& face_index_;
    vector<Vector2D> const& face_velocity_;
  };
}

double CalcTimeStep(Tessellation const& tessellation,
		    vector<Primitive> const& cells,
		    vector<Vector2D> const& facevelocity,
		    HydroBoundaryConditions const& hbc,
		    double time,
		    vector<CustomEvolution*> const& evolve)
{
  bool first_time=true;
  double dt=0;
  for(int i = 0;i<tessellation.GetPointNo(); i++)
    {
      if(irrelevant_for_time_step(evolve,i))
	continue;
      const vector<Vector2D> face_vel = initialize_vector<Vector2D>
	(FaceVelocityInitializer(tessellation.GetCellEdges(i),
				 facevelocity));
      const double dt_temp = calc_dt_temp
	(tessellation, cells, face_vel,
	 hbc, time, i);
      if(first_time)
	{
	  first_time=false;
	  dt=dt_temp;
	}
      else
	dt = min(dt_temp, dt);
    }
  return dt;
}

namespace {
  void update_conserved_extensive_error(int edge_index,
					int cell_number)
  {
    UniversalError eo("Error in UpdateConservedExtensive: cell and edge are not mutual neighbors");
    eo.AddEntry("edge number",edge_index);
    eo.AddEntry("cell number",cell_number);
    throw eo;
  }
}

void UpdateConservedExtensive
(Tessellation const& tessellation,
 vector<Conserved> const& fluxes,
 double dt, 
 vector<Conserved>& conserved_extensive,
 HydroBoundaryConditions const& boundaryconditions)
{
  for(int i=0;i<tessellation.GetPointNo();++i){
    const vector<int> cell_edge_indices = tessellation.GetCellEdges(i);
    if(!boundaryconditions.IsGhostCell(i,tessellation)){
      for(int j=0;j<(int)cell_edge_indices.size();++j){
	const int edge_index = cell_edge_indices[j];
	const Edge edge = tessellation.GetEdge(edge_index);
	const Conserved delta = dt*edge.GetLength()*fluxes[edge_index];
	if(i==edge.GetNeighbor(0))
	  conserved_extensive[i] -= delta;
	else if(i==edge.GetNeighbor(1))
	  conserved_extensive[i] += delta;
	else
	  update_conserved_extensive_error(edge_index,i);
      }
    }
  }
}

namespace {
  class NewPointPosition: public VectorTermInitializer<Vector2D>
  {
  public:

    NewPointPosition(Tessellation const& tess,
		     vector<Vector2D> const& point_velocity,
		     double dt):
      tess_(tess),
      point_velocity_(point_velocity),
      dt_(dt) {}

    Vector2D operator()(int n) const
    {
      return tess_.GetMeshPoint(n)+dt_*point_velocity_[n];
    }

    int getLength(void) const
    {
      return tess_.GetPointNo();
    }

  private:
    Tessellation const& tess_;
    vector<Vector2D> const& point_velocity_;
    const double dt_;
  };
}

void MoveMeshPoints(vector<Vector2D> const& pointvelocity,
		    double dt, Tessellation& tessellation)
{
  tessellation.Update(initialize_vector<Vector2D>
		      (NewPointPosition(tessellation,
					pointvelocity,
					dt)));
}

void UpdateConservedIntensive(Tessellation const& tessellation,
			      vector<Conserved> const& conservedextensive,
			      vector<Conserved>& conservedintensive)
{
  for(int i = 0; i<tessellation.GetPointNo(); i++){
    conservedintensive[i] = conservedextensive[i]/
      tessellation.GetVolume(i);
  }
}

namespace {

  Conserved density_floor_correction(double density_min,
				     EquationOfState const& eos,
				     Primitive const& old_cell)
  {
    const double mass = density_min;
    const double new_pressure = old_cell.Pressure*density_min/old_cell.Density;
    const double kinetic_energy = 0.5*pow(abs(old_cell.Velocity),2)*density_min;
    const double thermal_energy = density_min*
      eos.dp2e(density_min, new_pressure);
    const Vector2D momentum = density_min*old_cell.Velocity;
    return Conserved(mass,momentum,kinetic_energy+thermal_energy);    
  }

  Conserved calc_safe_conserved(Conserved const& raw,
				bool density_floor,
				double min_density,
				double min_pressure,
				Primitive const& old,
				EquationOfState const& eos)
  {
    Conserved res = raw;
    if(density_floor){
      if(res.Mass<min_density)
	res = density_floor_correction
	  (min_density, eos, old);
      const double kinetic_energy = 0.5*pow(abs(res.Momentum),2)/res.Mass;
      const double thermal_energy = res.Energy-kinetic_energy;
      const double pressure = eos.de2p(res.Mass,thermal_energy/res.Mass);
      if(pressure<min_pressure)
	res.Energy = kinetic_energy+res.Mass*eos.dp2e(res.Mass, min_pressure);
    }
    return res;
  }

  void update_primitives_rethrow(int cell_index,
				 UniversalError& eo)
  {
    eo.AddEntry("UpdatePrimitive data starts here",0);
    eo.AddEntry("cell index",(double)cell_index);
    throw eo;
  }

  Primitive regular_cell_evolve(Conserved const& intensive,
				bool density_floor,
				double min_density,
				double min_pressure,
				Primitive const& old,
				EquationOfState const& eos)
  {
    const Conserved temp = calc_safe_conserved
      (intensive,density_floor, min_density,
       min_pressure, old, eos);
    return Conserved2Primitive(temp, eos);
  }
}

void UpdatePrimitives
(vector<Conserved> const& conservedintensive,
 EquationOfState const& eos,vector<Primitive>& cells,
 vector<CustomEvolution*> const& CellsEvolve,vector<Primitive> &old_cells,
 bool densityfloor,double densitymin,double pressuremin,Tessellation const&
 tess,double time,vector<vector<double> > const& tracers)
{
  for(int i=0;i < tess.GetPointNo(); i++){
    try
      {
	if(CellsEvolve[i]==0)
	  cells[i] = regular_cell_evolve
	    (conservedintensive[i], densityfloor,
	     densitymin, pressuremin, old_cells[i], eos);
	else
	  cells[i]=CellsEvolve[i]->UpdatePrimitive
	    (conservedintensive,
	     eos,old_cells,i,tess,time,tracers);
      }
    catch(UniversalError& eo){
      update_primitives_rethrow(i,eo);
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
  const Primitive rotated_left = RotatePrimitive(normaldir, paraldir, left);
  const Primitive rotated_right = RotatePrimitive(normaldir, paraldir, right);
  const double normal_speed = Projection(edge_velocity,normaldir);
  const Conserved res = rs.Solve(rotated_left, rotated_right, normal_speed);
  return RotateFluxBack(res, normaldir, paraldir);
}

namespace {

  Conserved calc_single_flux_in_bulk(Tessellation const& tess,
				     Edge const& edge,
				     SpatialReconstruction const& interpolation,
				     vector<Primitive> const& cells,
				     Vector2D const& face_velocity,
				     RiemannSolver const& rs,
				     double dt)
  {
    const Vector2D normal_dir = 
      tess.GetMeshPoint(edge.GetNeighbor(1))-
      tess.GetMeshPoint(edge.GetNeighbor(0));

    const Vector2D paral_dir = 
      edge.GetVertex(1) - edge.GetVertex(0);

    const Primitive left = interpolation.Interpolate
      (tess,cells,dt,edge,0,InBulk,face_velocity);

    const Primitive right = interpolation.Interpolate
      (tess,cells,dt,edge,1,InBulk,face_velocity);

    return FluxInBulk(normal_dir, paral_dir,
		      left, right,
		      face_velocity, rs);
		      
  }

  int choose_special_cell_index
  (vector<CustomEvolution*> const& ce_list,
   int n0, int n1)
  {
    if(ce_list[n0]&&ce_list[n1])
      throw UniversalError("Error in choose_special_cell_index: Both sides are special cells");
    else if(!ce_list[n0]&&!ce_list[n1])
      throw UniversalError("Error in choose_special_cell_index: Both sides are regular cells");
    else if(ce_list[n0]&&!ce_list[n1])
      return n0;
    else if(!ce_list[n0]&&ce_list[n1])
      return n1;
    else
      throw UniversalError("Error in choose_special_cell_index: Something has gone terribly wrong if you've gotten here");
  }
 
  void calc_fluxes_rethrow(UniversalError& eo,
			   int edge_index,
			   Tessellation const& tess)
  {
    eo.AddEntry("Error in CalcFlux",0);
    eo.AddEntry("edge index",edge_index);
    const Edge edge = tess.GetEdge(edge_index);
    eo.AddEntry("edge x1 location",edge.get_x(0));
    eo.AddEntry("edge y1 location",edge.get_y(0));
    eo.AddEntry("edge x2 location",edge.get_x(1));
    eo.AddEntry("edge y2 location",edge.get_y(1));
    throw eo;
  }

  bool should_skip_flux_calculation(CustomEvolution* c1,
				    CustomEvolution* c2)
  {
    if(c1&&c2){
      if(c1->flux_indifferent()&&c2->flux_indifferent())
	return true;
    }
    return false;
  }
}

void CalcFluxes
(Tessellation const& tessellation,
 vector<Primitive> const& cells,
 double dt,
 double time,
 SpatialReconstruction& interpolation,
 vector<Vector2D> const& facevelocity,
 HydroBoundaryConditions const& boundaryconditions,
 RiemannSolver const& rs,
 vector<Conserved>& fluxes,
 vector<CustomEvolution*> const& CellsEvolve,
 vector<vector<double> > const& tracers)
{
  fluxes.resize(tessellation.GetTotalSidesNumber());

  try
    {
      interpolation.Prepare(tessellation,cells,tracers,dt,time);
    }
  catch(UniversalError& eo)
    {
      eo.AddEntry("Error in CalcFlux Interpolation prepare",0);
      throw;
    }
  for(int i = 0; i < tessellation.GetTotalSidesNumber(); i++)
    {
      try
	{
	  const Edge edge = tessellation.GetEdge(i);
	  const int n1=edge.GetNeighbor(1);
	  const int n0=edge.GetNeighbor(0);
	  if(!boundaryconditions.IsBoundary(edge,tessellation))
	    {
	      if(!CellsEvolve[n0]&&!CellsEvolve[n1])
		fluxes[i] = calc_single_flux_in_bulk
		  (tessellation, edge, interpolation,
		   cells,facevelocity[i],rs,dt);
	      else
		{
		  if(!should_skip_flux_calculation(CellsEvolve[n0],
						   CellsEvolve[n1])){
		    const int ns = choose_special_cell_index
		      (CellsEvolve, n0, n1);
		    fluxes.at(i) = CellsEvolve.at(ns)->CalcFlux
		      (tessellation, cells,
		       dt,interpolation,edge,facevelocity[i],
		       rs,ns,boundaryconditions,time,tracers);
		  }
		}
	    }
	  else
	    {
	      // Boundaries
	      fluxes[i] = boundaryconditions.CalcFlux
		(tessellation, cells, facevelocity[i], edge,
		 interpolation,dt,time);
	    }
	}
      catch(UniversalError& eo)
	{
	  calc_fluxes_rethrow(eo,i,tessellation);
	}
    }
}

void ExternalForceContribution
(Tessellation const& tess,
 vector<Primitive> const& cells,
 SourceTerm& force,
 double t,
 double dt,
 vector<Conserved>& conserved_extensive,
 HydroBoundaryConditions const& hbc,vector<Conserved> const& fluxes,
 vector<Vector2D> const& point_velocity,vector<double> &g,bool coldflows_flag,
 vector<vector<double> > &tracers_extensive)
{
  vector<double> dtracer;
  if(!tracers_extensive.empty())
    dtracer.assign(tracers_extensive[0].size(),0);
  if(coldflows_flag)
    g.resize(point_velocity.size());
  for(int i=0;i<tess.GetPointNo();++i)
    {
      if(!hbc.IsGhostCell(i,tess))
	{
	  const Conserved cons(force.Calculate
			       (tess,cells,i,
				fluxes,point_velocity,hbc,tracers_extensive,
				dtracer,t,dt));
	  conserved_extensive[i]+=dt*cons;
	  if(!tracers_extensive.empty())
	    {
	      tracers_extensive[i]=tracers_extensive[i]+dt*dtracer;
	    }
	  if(coldflows_flag)
	    g[i]=abs(cons.Momentum)/(cells[i].Density*tess.GetVolume(i));
	}
    }
}


vector<Vector2D> calc_point_velocities
(Tessellation const& tess,
 vector<Primitive> const& cells,
 PointMotion& point_motion,
 double time)
{
  return point_motion.calcAllVelocities(tess,cells,time);
}

namespace {
  vector<Conserved> calc_fluxes
  (Tessellation const& tess,
   vector<Primitive> const& cells,
   double dt,
   double time,
   SpatialReconstruction& interpolation,
   vector<Vector2D> const& edge_velocities,
   HydroBoundaryConditions const& hbc,
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
(Tessellation const& tess)
{
  vector<Vector2D> res(tess.GetPointNo());
  for(int i=0;i<(int)tess.GetPointNo();++i)
    res[i] = tess.GetMeshPoint(i);
  return res;
}

double TimeAdvance2mid
(Tessellation& tess,vector<Primitive> &cells,
 PointMotion& point_motion,HydroBoundaryConditions const& hbc,
 SpatialReconstruction& interpolation,RiemannSolver const& rs,
 EquationOfState const& eos,SourceTerm& force,double time,double cfl,
 double endtime,vector<CustomEvolution*> const& CellsEvolve,
 vector<vector<double> > &tracers,
 double dt_external,bool traceflag,bool coldflows_flag,
 double as,double bs,bool densityfloor,double densitymin,
 double pressuremin,bool EntropyCalc)
{
  boost::scoped_ptr<Tessellation> sp_tess(tess.clone());
  Tessellation& tess_temp = *sp_tess.get();

  //do half time step
  vector<Vector2D> point_velocities = calc_point_velocities
    (tess_temp, cells, point_motion, time);

  vector<Vector2D> edge_velocities =tess_temp.calc_edge_velocities
    (&hbc,
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
      const int n=tess.GetPointNo();
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

  vector<Primitive> new_cells(tess_temp.GetPointNo());
  UpdatePrimitives(intensive, eos, new_cells,CellsEvolve,cells,densityfloor,
		   densitymin,pressuremin,tess_temp,time+0.5*dt,tracers);
  if(traceflag)
    {
      MakeTracerIntensive(tracers,tracer_extensive,tess_temp,new_cells);
    }

  // End half step

  point_velocities = calc_point_velocities(tess_temp,new_cells,
					   point_motion, time+0.5*dt);
  edge_velocities = tess_temp.calc_edge_velocities(&hbc,point_velocities,time);

  fluxes = calc_fluxes
    (tess_temp, new_cells, dt, time, interpolation,
     edge_velocities, hbc, rs,CellsEvolve,tracers);

  if(coldflows_flag&&EntropyCalc)
    {
      const int n=tess.GetPointNo();
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

namespace {
  pair<CustomEvolution*,int> choose_time_evolve_method
  (pair<CustomEvolution*, int> const& opt_1,
   pair<CustomEvolution*, int> const& opt_2)
  {
    if(opt_1.first&&
       opt_2.first&&
       !(opt_1.first->flux_indifferent()&&
	 opt_2.first->flux_indifferent())){
      UniversalError eo("Bad inner boundary, edge doesn't have same cellevolve conditions on both sides");
      eo.AddEntry("first cell",opt_1.second);
      eo.AddEntry("second cell",opt_2.second);
      throw eo;
    }
    else if(opt_1.first&&!opt_2.first)
      return opt_1;
    else if(opt_2.first&&!opt_1.first)
      return opt_2;
    else
      return pair<CustomEvolution*,int>((CustomEvolution*)0,0);
  }

  class ConditionalPlusMinus: public BinaryOperation<double>
  {
  public:

    ConditionalPlusMinus(bool flag):
      flag_(flag) {}

    double operator()(double const& x, double const& y) const
    {
      if(flag_)
	return x+y;
      else
	return x-y;
    }

  private:
    const bool flag_;
  };

  class ScalarMultiply: public UnaryOperation<double>
  {
  public:

    ScalarMultiply(double scalar):
      scalar_(scalar) {}

    double operator()(double const& x) const
    {
      return x*scalar_;
    }

  private:
    const double scalar_;
  };
}

vector<vector<double> > CalcTraceChange
(vector<vector<double> > const& old_trace,vector<Primitive> const& cells,
 Tessellation const& tess,vector<Conserved> const& fluxes,double dt,
 HydroBoundaryConditions const& hbc, SpatialReconstruction const& interp,
 double time,vector<CustomEvolution*> const& CellsEvolve,
 vector<Vector2D> const& edge_velocities)
{
  if(old_trace.empty())
    return vector<vector<double> >();
  const int n = (int)old_trace.size();
  const int dim = (int)old_trace[0].size();
  vector<vector<double> > res(n,vector<double>(dim,0));
  for(int i=0;i<n;++i)
    {
      const vector<int> cell_edge_indices=tess.GetCellEdges(i);
      const int Nedges=int(cell_edge_indices.size());
      for(int k=0;k<Nedges;++k)
	{
	  const int edge_index = cell_edge_indices[k];
	  const Edge edge = tess.GetEdge(edge_index);
	  const double dm=fluxes[edge_index].Mass;
	  const int n1=edge.GetNeighbor(1);
	  const int n0=edge.GetNeighbor(0);
	  if(!hbc.IsBoundary(edge,tess))
	    {
	      const pair<CustomEvolution*,int> evol_method =
		choose_time_evolve_method
		(pair<CustomEvolution*,int>(CellsEvolve[n0],n0),
		 pair<CustomEvolution*,int>(CellsEvolve[n1],n1));
	      if(evol_method.first){
		const vector<double> temp = evol_method.first->CalcTracerFlux
		  (tess,cells,old_trace,dm,edge,evol_method.second,dt,
		   time,interp,edge_velocities[edge_index]);
		transform(res[i].begin(),
			  res[i].end(),
			  temp.begin(),
			  res[i].begin(),
			  plus<double>());
		continue;
	      }
	    }

	  if(hbc.IsBoundary(edge,tess))
	    {
	      const vector<double> temp=hbc.CalcTracerFlux
		(tess,cells,old_trace,dm,edge,i,dt,time,interp,
		 edge_velocities[edge_index]);
	      res[i] = binary_unite(res[i],temp,
				 ConditionalPlusMinus(i==n1));
	    }
	  else
	    {
	      const vector<double> temp2=interp.interpolateTracers
	      (tess,cells,old_trace,dt,edge,dm<0,InBulk,
	       edge_velocities[edge_index]);
	      const vector<double> temp = apply_to_each_term
	      (temp2,ScalarMultiply(dm*dt*edge.GetLength()));
	      res[i] = binary_unite(res[i],temp,
				    ConditionalPlusMinus(i==n1));
	    }
	}
    }
	    
  return res;
}

vector<double> GetMaxKineticEnergy(Tessellation const& tess,vector<Primitive> const&
				   cells,vector<CustomEvolution*> const& customevolve)
{
  const int n=tess.GetPointNo();
  vector<double> res;
  res.resize(n);
  double e;
  for(int j=0;j<n;++j)
    {
      e=0;
      if(!NearBoundary(j,tess,customevolve))
	{	
	  vector<int> neigh=tess.GetNeighbors(j);
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

vector<double> GetForceEnergy(Tessellation const& tess,
			      vector<double> const& g)
{
  vector<double> res;
  int n=int(g.size());
  res.resize(n);
  for(int i=0;i<n;++i)
    res[i]=g[i]*tess.GetWidth(i);
  return res;
}

void FixPressure(vector<Conserved> &intensive,vector<vector<double> > const& entropy,
		 EquationOfState const& eos,vector<double> const& Ek,
		 vector<double> const& Ef,double as,double bs,vector<CustomEvolution*>
		 const& customevolve,Tessellation const& tess,vector<Conserved> &extensive,
		 vector<Primitive> const& cells,HydroBoundaryConditions const& hbc,
		 double time)
{
  int n=int(intensive.size());
  double Et,Ek2;
  double temp;
  for(int i=0;i<n;++i){
    if(customevolve[i]==0){
      //Make intensive
      temp=entropy[i][0]/(tess.GetVolume(i)*intensive[i].Mass);
      Ek2=0.5*pow(abs(intensive[i].Momentum)/intensive[i].Mass,2);
      Et=intensive[i].Energy/intensive[i].Mass-Ek2;
      if((Et<as*Ek[i])||(Et<bs*Ef[i])){
	if(~IsShockedCell(tess,i,cells,hbc,time)||Et<0){
	  Et=eos.dp2e(intensive[i].Mass,
		      eos.sd2p(temp,intensive[i].Mass));
	  intensive[i].Energy=intensive[i].Mass*(Et+Ek2);
	  extensive[i].Energy=tess.GetVolume(i)*intensive[i].Energy;
	}
      }
    }
  }
}

bool NearBoundary(int index,Tessellation const& tess,
		  vector<CustomEvolution*> const& /*customevolve*/)
{
  vector<int> neigh=tess.GetNeighbors(index);
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
			 Tessellation const& tess,vector<Primitive> const& cells,vector<vector<double> >
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
      mass=tess.GetVolume(i)*cells[i].Density;
      result[i].resize(dim);
      for(int j=0;j<dim;++j)
	result[i][j]=tracer[i][j]*mass;
    }
  return;
}

void MakeTracerIntensive(vector<vector<double> > &tracer,vector<vector<double> >
			 const& tracer_extensive,Tessellation const& tess,vector<Primitive> const& cells)
{
  int n=int(tracer.size());
  if(n==0)
    return;
  int dim=int(tracer[0].size());
  double mass;
  for(int i=0;i<n;++i)
    {
      mass=tess.GetVolume(i)*cells[i].Density;
      for(int j=0;j<dim;++j)
	tracer[i][j]=tracer_extensive[i][j]/mass;
    }
  return;
}

void UpdateTracerExtensive(vector<vector<double> > &tracerextensive,
			   vector<vector<double> > const& tracerchange,vector<CustomEvolution*> const&
			   CellsEvolve,vector<Primitive> const& cells,Tessellation const& tess,
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

void TracerResetCalc
(double alpha,SpatialDistribution const& originalD,
 SpatialDistribution const& originalP,SpatialDistribution const& originalVx,
 SpatialDistribution const& originalVy, vector<Primitive> &cells,
 Tessellation const& tess,vector<vector<double> > &tracer,
 int tracerindex,EquationOfState const& eos,int InnerNum)
{
  const int n = tess.GetPointNo();
  if(n<1)
    return;
  Vector2D velocity;
  if(tracer.empty())
    {
      for(int i=0;i<InnerNum;++i)
	{
	  velocity.Set(originalVx.EvalAt(tess.GetCellCM(i)),
		       originalVy.EvalAt(tess.GetCellCM(i)));
	  cells[i]=CalcPrimitive(originalD.EvalAt(tess.GetCellCM(i)),
				 originalP.EvalAt(tess.GetCellCM(i)),velocity,eos);
	}
      return;	
    }
  if(tracerindex>=(int)tracer[0].size()||tracerindex<0)
    throw UniversalError("Error in tracerReset, wrong dimension for tracer");
  for(int i=0;i<n;++i)
    {
      if((tracer[i][tracerindex]<alpha)||(i<InnerNum))//*cells[i].Density)
	{
	  velocity.Set(originalVx.EvalAt(tess.GetCellCM(i)),
		       originalVy.EvalAt(tess.GetCellCM(i)));
	  cells[i]=CalcPrimitive(originalD.EvalAt(tess.GetCellCM(i)),
				 originalP.EvalAt(tess.GetCellCM(i)),velocity,eos);
	  if(tracer[i][tracerindex]<0)
	    tracer[i][tracerindex]=0;
	  if(i<InnerNum)
	    tracer[i][tracerindex]=0;
	}
    }
  return;
}

void GetPointToRemove(Tessellation const& tess,Vector2D const& point,
		      double R,vector<int> & PointToRemove,int Inner)
{
  int n=tess.GetPointNo();
  PointToRemove.clear();
  bool test;
  for(int i=Inner;i<n;++i)
    {
      // Check if point is completly engulfed
      test=true;
      vector<int> neigh=tess.GetNeighbors(i);
      for(int j=0;j<(int)neigh.size();++j)
	if(neigh[j]>=Inner)
	  test=false;
      // Is point inside a radius?
      if(abs(point-tess.GetMeshPoint(i))<R||test)
	PointToRemove.push_back(i);
    }
  return;
}

namespace {
  Vector2D GetReflectedPoint(Tessellation const& tess,int point,
			     Edge const& edge)
  {
    Vector2D MeshPoint=tess.GetMeshPoint(point);
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

bool IsShockedCell(Tessellation const& tess,int index,
		   vector<Primitive> const& cells,
		   HydroBoundaryConditions const& hbc,
		   double time)
{
  double DivV=0;
  vector<int> edges_loc=tess.GetCellEdges(index);
  vector<Edge> edges(edges_loc.size());
  for(int i=0;i<(int)edges.size();++i)
    edges[i]=tess.GetEdge(edges_loc[i]);

  // Calculate gradiant by gauss theorem
  double vx_i = cells[index].Velocity.x;
  double vy_i = cells[index].Velocity.y;
  int other;
  Vector2D center=tess.GetMeshPoint(index);
  Vector2D c_ij,r_ij;
  for(int i=0;i<(int)edges.size();i++)
    {
      if(hbc.IsBoundary(edges[i],tess))
	{
	  Primitive other2=hbc.GetBoundaryPrimitive(edges[i],tess,cells,
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
						 tess.GetMeshPoint(other_point)+center);
	      r_ij = center-tess.GetMeshPoint(other_point);
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
	    0.5*(tess.GetMeshPoint(other)+center);
	  r_ij = center-tess.GetMeshPoint(other);
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
