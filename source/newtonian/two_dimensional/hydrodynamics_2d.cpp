#include "hydrodynamics_2d.hpp"
#include "../../tessellation/calc_face_vertex_velocity.hpp"

using std::max;

FaceVertexVelocityCalculator::FaceVertexVelocityCalculator
(const Tessellation& tess,
 const vector<Vector2D>& point_velocities,
 const Vector2D std::pair<Vector2D,Vector2D>::* const member,
 const vector<Vector2D>& control,
 const HydroBoundaryConditions& hbc):
  tess_(tess),
  point_velocities_(point_velocities),
  member_(member),
  control_(control),
  hbc_(hbc) {}

size_t FaceVertexVelocityCalculator::getLength(void) const
{
  return static_cast<size_t>(tess_.GetTotalSidesNumber());
}

Vector2D FaceVertexVelocityCalculator::operator()(size_t i) const
{
  const Edge& edge = tess_.GetEdge(static_cast<int>(i));
  if(hbc_.IsBoundary(edge,tess_))
    return control_[i];
  return calc_face_vertex_velocity
    (tess_.GetMeshPoint(edge.neighbors.first),
     point_velocities_[static_cast<size_t>(edge.neighbors.first)],
     tess_.GetMeshPoint(edge.neighbors.second),
     point_velocities_[static_cast<size_t>(edge.neighbors.second)],
     edge.vertices.*member_);
}

int get_other_index(const Edge& edge, const int index)
{
  if(edge.neighbors.first==index && edge.neighbors.second!=index)
    return edge.neighbors.second;
  else if(edge.neighbors.second==index && edge.neighbors.first!=index)
    return edge.neighbors.first;
  else
    throw UniversalError("Something went wrong in Hydrodynamics2D::get_other_index");
}

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
    return CalcPrimitive(density(r),
			 pressure(r),
			 Vector2D(xvelocity(r),
				  yvelocity(r)),
			 eos);
  }

  class CellInitializer: public LazyList<Primitive>
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

    Primitive operator[](size_t n) const
    {
      return initialize_single_cell(tess_,
				    static_cast<int>(n),
				    cm_flag_,
				    density_,
				    pressure_,
				    xvelocity_,
				    yvelocity_,
				    eos_);
    }

    size_t size(void) const
    {
      return static_cast<size_t>(tess_.GetPointNo());
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
  return serial_generate
    (CellInitializer
     (tess, cm_value, density,
      pressure, xvelocity, yvelocity, eos));
}

namespace {
  class IntensiveInitializer: public LazyList<Conserved>
  {
  public:

    IntensiveInitializer(vector<Primitive> const& cells):
      cells_(cells) {}

    Conserved operator[](size_t n) const
    {
      return Primitive2Conserved(cells_[n]);
    }

    size_t size(void) const
    {
      return cells_.size();
    }

  private:
    vector<Primitive> const& cells_;
  };
}

vector<Conserved> CalcConservedIntensive
(vector<Primitive> const& cells)
{
  return serial_generate(IntensiveInitializer(cells));
}

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

  class ExtensiveInitializer: public LazyList<Conserved>
  {
  public:

    ExtensiveInitializer(const vector<Conserved>& intensive,
			 const Tessellation& tess,
			 const PhysicalGeometry& pg):
      intensive_(intensive), tess_(tess), pg_(pg)  {}

    Conserved operator[](size_t n) const
    {
      return pg_.calcVolume(serial_generate(CellEdgesGetter(tess_,static_cast<int>(n))))*
	intensive_[n];
    }

    size_t size(void) const
    {
      return static_cast<size_t>(tess_.GetPointNo());
    }

  private:
    const vector<Conserved>& intensive_;
    const Tessellation& tess_;
    const PhysicalGeometry& pg_;
  };
}

vector<Conserved> CalcConservedExtensive
(const vector<Conserved>& cons_int,
 const Tessellation& tess,
 const PhysicalGeometry& pg)
{
  return serial_generate(ExtensiveInitializer(cons_int, tess, pg));
}

double TimeStepForCell(Primitive const& cell,double width,
		       vector<Vector2D> const& face_velocites)
{
  double max_fv=0;
  for(size_t i=0;i<face_velocites.size();++i)
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
      return max(abs(face_velocities[static_cast<size_t>(i)]-cell.Velocity),
		 abs(face_velocities[static_cast<size_t>(i)]-
		     hbc.GetBoundaryPrimitive(edge,
					      tess,
					      cells,
					      time).Velocity));
    else
      return abs(face_velocities[static_cast<size_t>(i)]-cell.Velocity);
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
  for(int i=0;i<static_cast<int>(face_velocities.size());++i)
    {
      const Edge& edge=tess.GetEdge(edge_index[static_cast<size_t>(i)]);
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
      if(cevolve[static_cast<size_t>(i)]!=0)
	if(!cevolve[static_cast<size_t>(i)]->TimeStepRelevant())
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
      return TimeStepForCell(cells[static_cast<size_t>(i)],tess.GetWidth(i),face_vel);
    else
      return TimeStepForCellBoundary(cells[static_cast<size_t>(i)],
				     cells,
				     tess.GetWidth(i),
				     face_vel,
				     tess,
				     hbc,
				     i,
				     time);
  }

  class FaceVelocityInitializer: public LazyList<Vector2D>
  {
  public:

    FaceVelocityInitializer(vector<int> const& face_index,
			    vector<Vector2D> const& face_velocity):
      face_index_(face_index),
      face_velocity_(face_velocity) {}

    Vector2D operator[](size_t n) const
    {
      return face_velocity_[static_cast<size_t>(face_index_[n])];
    }

    size_t size(void) const
    {
      return face_index_.size();
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
      const vector<Vector2D> face_vel = serial_generate
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
#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
  __attribute__((noreturn))
#endif
  void update_conserved_extensive_error
  (int edge_index, int cell_number)
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
 HydroBoundaryConditions const& boundaryconditions,
 vector<double> const& lengthes)
{
  for(int i=0;i<tessellation.GetPointNo();++i){
    const vector<int> cell_edge_indices = tessellation.GetCellEdges(i);
    if(!boundaryconditions.IsGhostCell(i,tessellation)){
      for(int j=0;j<static_cast<int>(cell_edge_indices.size());++j){
	const int edge_index = cell_edge_indices[static_cast<size_t>(j)];
	const Edge& edge = tessellation.GetEdge(edge_index);
	const Conserved delta = dt*lengthes[static_cast<size_t>(edge_index)]*fluxes[static_cast<size_t>(edge_index)];
	if(i==edge.neighbors.first)
	  conserved_extensive[static_cast<size_t>(i)] -= delta;
	else if(i==edge.neighbors.second)
	  conserved_extensive[static_cast<size_t>(i)] += delta;
	else
	  update_conserved_extensive_error(edge_index,i);
      }
    }
  }
}

namespace {
  class NewPointPosition: public LazyList<Vector2D>
  {
  public:

    NewPointPosition(Tessellation const& tess,
		     vector<Vector2D> const& point_velocity,
		     double dt):
      tess_(tess),
      point_velocity_(point_velocity),
      dt_(dt) {}

    Vector2D operator[](size_t n) const
    {
      return tess_.GetMeshPoint(static_cast<int>(n))+dt_*point_velocity_[n];
    }

    size_t size(void) const
    {
      return static_cast<size_t>(tess_.GetPointNo());
    }

  private:
    Tessellation const& tess_;
    vector<Vector2D> const& point_velocity_;
    const double dt_;
  };
}

void MoveMeshPoints(vector<Vector2D> const& pointvelocity,
		    double dt, Tessellation& tessellation,vector<Vector2D> oldpoints)
{
  if(oldpoints.empty())
    {
      tessellation.Update(serial_generate
			  (NewPointPosition(tessellation,
					    pointvelocity,
					    dt)));
    }
  else
    {
      for(size_t i=0;i<oldpoints.size();++i)
	oldpoints[i]+=pointvelocity[i]*dt;
      tessellation.Update(oldpoints);
    }
}

#ifdef RICH_MPI
void MoveMeshPoints(vector<Vector2D> const& pointvelocity,
		    double dt, Tessellation& tessellation,

		    Tessellation const& vproc,
		    vector<Vector2D> oldpoints)
{
  if(oldpoints.empty())
    tessellation.Update(serial_generate
			(NewPointPosition(tessellation,
					  pointvelocity,
					  dt)),vproc);
  else
    {
      for(size_t i=0;i<oldpoints.size();++i)
	oldpoints[i]+=pointvelocity[i]*dt;
      tessellation.Update(oldpoints,vproc);
    }
}
#endif // RICH_MPI

namespace {
  class IntensiveCalculator: public LazyList<Conserved>
  {
  public:

    IntensiveCalculator(const Tessellation& tess,
			const vector<Conserved>& extensive,
			const PhysicalGeometry& pg):
      tess_(tess), extensive_(extensive), pg_(pg) {}

    size_t size(void) const
    {
      return extensive_.size();
    }

    Conserved operator[](size_t i) const
    {
      return extensive_[i]/
	pg_.calcVolume(serial_generate(CellEdgesGetter(tess_,static_cast<int>(i))));
    }

  private:
    const Tessellation& tess_;
    const vector<Conserved>& extensive_;
    const PhysicalGeometry& pg_;
  };
}

vector<Conserved> calc_conserved_intensive
(const Tessellation& tess,
 const vector<Conserved>& extensive,
 const PhysicalGeometry& pg)
{
  return serial_generate(IntensiveCalculator(tess,extensive,pg));
}

void UpdateConservedIntensive(Tessellation const& tessellation,
			      vector<Conserved> const& conservedextensive,
			      vector<Conserved>& conservedintensive)
{
  conservedintensive.resize(conservedextensive.size());
  for(int i = 0; i<tessellation.GetPointNo(); i++){
    conservedintensive[static_cast<size_t>(i)] = conservedextensive[static_cast<size_t>(i)]/
      tessellation.GetVolume(i);
  }
}

namespace {

  std::pair<Conserved,bool> calc_safe_conserved(Conserved const& raw,
						bool density_floor,
						double min_density,
						double min_pressure,
						Primitive const& old,
						EquationOfState const& eos)
  {
    std::pair<Conserved, bool> res;
    res.first = raw;
    if(density_floor)
      {
	if (res.first.Mass < min_density)
	  {
	    res.first.Mass = min_density;
	    res.first.Momentum = old.Velocity*min_density;
	    const double kinetic_energy = 0.5*pow(abs(res.first.Momentum / res.first.Mass), 2);
	    res.first.Energy = res.first.Mass*kinetic_energy
	      + res.first.Mass*eos.dp2e(res.first.Mass, min_pressure);
	    res.second = true;
	  }
	const double kinetic_energy = 0.5*pow(abs(res.first.Momentum / res.first.Mass), 2);
	const double thermal_energy = res.first.Energy / res.first.Mass - kinetic_energy;
	const double pressure = eos.de2p(res.first.Mass, thermal_energy);
	if (pressure < min_pressure || res.second)
	  {
	    res.first.Energy = res.first.Mass*kinetic_energy 
	      + res.first.Mass*eos.dp2e(res.first.Mass, min_pressure);
	    res.second = true;
	  }
      }
    return res;
  }

#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
  __attribute__((noreturn))
#endif
  void update_primitives_rethrow(int cell_index,
				 UniversalError& eo)
  {
    eo.AddEntry("UpdatePrimitive data starts here",0);
    eo.AddEntry("cell index",static_cast<double>(cell_index));
    throw eo;
  }

  std::pair<Primitive,bool> regular_cell_evolve(Conserved const& intensive,
						bool density_floor,
						double min_density,
						double min_pressure,
						Primitive const& old,
						EquationOfState const& eos)
  {
    const std::pair<Conserved, bool> temp = calc_safe_conserved
      (intensive,density_floor, min_density,
       min_pressure, old, eos);
    return std::pair<Primitive,bool> (Conserved2Primitive(temp.first, eos),temp.second);
  }
}

vector<bool> UpdatePrimitives
(vector<Conserved> const& conservedintensive,
 EquationOfState const& eos,vector<Primitive>& cells,
 vector<CustomEvolution*> const& CellsEvolve,vector<Primitive> &old_cells,
 bool densityfloor,double densitymin,double pressuremin,Tessellation const&
 tess,double time,vector<vector<double> > const& tracers)
{
  cells.resize(conservedintensive.size());
  vector<bool> bres(cells.size(), false);
  for(int i=0;i < tess.GetPointNo(); i++)
    {
      try
	{
	  if (CellsEvolve[static_cast<size_t>(i)] == 0 || CellsEvolve[static_cast<size_t>(i)]->DensityFloorRelevant())
	    {
	      Primitive old_cell = densityfloor ? old_cells[static_cast<size_t>(i)] : cells[static_cast<size_t>(i)];
	      std::pair<Primitive, bool > res = regular_cell_evolve
		(conservedintensive[static_cast<size_t>(i)], densityfloor,
		 densitymin, pressuremin, old_cell, eos);
	      if (CellsEvolve[static_cast<size_t>(i)] != 0 && !res.second)
		cells[static_cast<size_t>(i)] = CellsEvolve[static_cast<size_t>(i)]->UpdatePrimitive
		  (conservedintensive,eos, old_cells, i, tess, time, tracers);
	      else
		cells[static_cast<size_t>(i)] = res.first;
	      bres[static_cast<size_t>(i)] = res.second;
	    }
	  else
	    cells[static_cast<size_t>(i)]=CellsEvolve[static_cast<size_t>(i)]->UpdatePrimitive
	      (conservedintensive,
	       eos,old_cells,i,tess,time,tracers);
	}
      catch(UniversalError& eo)
	{
	  eo.AddEntry("x momentum per unit volume",conservedintensive[static_cast<size_t>(i)].Momentum.x);
	  eo.AddEntry("y momentum per unit volume",conservedintensive[static_cast<size_t>(i)].Momentum.y);
	  eo.AddEntry("thermal energy per unit mass",conservedintensive[static_cast<size_t>(i)].Energy);
	  eo.AddEntry("Cell volume",tess.GetVolume(i));
	  eo.AddEntry("Cell x location",tess.GetMeshPoint(i).x);
	  eo.AddEntry("Cell y location",tess.GetMeshPoint(i).y);
#ifdef RICH_MPI
	  eo.AddEntry("Error in CPU",static_cast<double>(get_mpi_rank()));
#endif
	  update_primitives_rethrow(i,eo);
	}
    }
  return bres;
}

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
  const Conserved res = rs(rotated_left, rotated_right, normal_speed);
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
      tess.GetMeshPoint(edge.neighbors.second)-
      tess.GetMeshPoint(edge.neighbors.first);

    const Vector2D paral_dir =
      edge.vertices.second - edge.vertices.first;

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
   const CustomEvolutionManager& cem,
   int n0, int n1)
  {
    if(ce_list[static_cast<size_t>(n0)]&&ce_list[static_cast<size_t>(n1)]){
      if(ce_list[static_cast<size_t>(n0)]==ce_list[static_cast<size_t>(n1)])
	return n0;
      const size_t priority_0 = cem.getIndex(ce_list[static_cast<size_t>(n0)]);
      const size_t priority_1 = cem.getIndex(ce_list[static_cast<size_t>(n1)]);
      assert(priority_0 != priority_1 &&
	     "Two methods have the same priority");
      return priority_0 > priority_1 ?
	n0 : n1;
    }
    else if(!ce_list[static_cast<size_t>(n0)]&&!ce_list[static_cast<size_t>(n1)])
      throw UniversalError("Error in choose_special_cell_index: Both sides are regular cells");
    else if(ce_list[static_cast<size_t>(n0)]&&!ce_list[static_cast<size_t>(n1)])
      return n0;
    else if(!ce_list[static_cast<size_t>(n0)]&&ce_list[static_cast<size_t>(n1)])
      return n1;
    else
      throw UniversalError("Error in choose_special_cell_index: Something has gone terribly wrong if you've gotten here");
  }

  /*
    void calc_fluxes_rethrow(UniversalError& eo,
    int edge_index,
    Tessellation const& tess)
    {
    eo.AddEntry("Error in CalcFlux",0);
    eo.AddEntry("edge index",edge_index);
    const Edge& edge = tess.GetEdge(edge_index);
    eo.AddEntry("edge x1 location",edge.vertices.first.x);
    eo.AddEntry("edge y1 location",edge.vertices.first.y);
    eo.AddEntry("edge x2 location",edge.vertices.second.x);
    eo.AddEntry("edge y2 location",edge.vertices.second.y);
    throw eo;
    }
  */
}

namespace {
  class InterpolationRelevancy: public LazyList<bool>
  {
  public:

    InterpolationRelevancy(const vector<CustomEvolution*>& ce):
      ce_(ce) {}

    size_t size(void) const
    {
      return ce_.size();
    }

    bool operator[](size_t i) const
    {
      if(ce_[i])
	return ce_[i]->isRelevantToInterpolation();
      else
	return true;
    }

  private:
    const vector<CustomEvolution*>& ce_;
  };

  class FluxCalculator: public LazyList<Conserved>
  {
  public:

    FluxCalculator(SpatialReconstruction& interp,
		   const Tessellation& tess,
		   const vector<Primitive>& cells,
		   const vector<vector<double> >& tracers,
		   const vector<CustomEvolution*>& ce,
		   const CustomEvolutionManager& cem,
		   double time, double dt,
		   const HydroBoundaryConditions& hbc,
		   const vector<Vector2D>& fv,
		   const RiemannSolver& rs):
      interp_(interp),
      tess_(tess),
      cells_(cells),
      tracers_(tracers),
      hbc_(hbc),
      ce_(ce),
      cem_(cem),
      fv_(fv),
      rs_(rs),
      time_(time),
      dt_(dt)
    {
      interp_.Prepare(tess,cells,tracers,
		      serial_generate(InterpolationRelevancy(ce)),
		      dt,time);
#ifndef RICH_MPI
      PeriodicGradExchange(interp_.GetGradients(),
			   tess_.GetDuplicatedPoints(),tess_.GetTotalPointNumber());
#else
      SendRecvGrad(interp_.GetGradients(),tess_.GetDuplicatedPoints(),
		   tess_.GetDuplicatedProcs(),tess_.GetGhostIndeces(),
		   tess_.GetTotalPointNumber());
#endif
    }

    size_t size(void) const
    {
      return static_cast<size_t>(tess_.GetTotalSidesNumber());
    }

    Conserved operator[](size_t i) const
    {
      const Edge& edge = tess_.GetEdge(static_cast<int>(i));
      const int n0 = edge.neighbors.first;
      const int n1 = edge.neighbors.second;
      if(!hbc_.IsBoundary(edge,tess_)){
	if(!ce_[static_cast<size_t>(n0)]&&!ce_[static_cast<size_t>(n1)])
	  return calc_single_flux_in_bulk(tess_,edge,interp_,cells_,fv_[i],rs_,dt_);
	else{
	  const int ns = choose_special_cell_index(ce_,cem_,n0,n1);
	  return ce_[static_cast<size_t>(ns)]->CalcFlux(tess_,cells_,dt_,interp_,edge,fv_[i],rs_,ns,hbc_,time_,tracers_);
	}	
      }
      else
	return hbc_.CalcFlux(tess_,cells_,fv_[i],edge,interp_,dt_,time_);
    }

  private:
    SpatialReconstruction& interp_;
    const Tessellation& tess_;
    const vector<Primitive>& cells_;
    const vector<vector<double> >& tracers_;
    const HydroBoundaryConditions& hbc_;
    const vector<CustomEvolution*>& ce_;
    const CustomEvolutionManager& cem_;
    const vector<Vector2D>& fv_;
    const RiemannSolver& rs_;
    const double time_;
    const double dt_;
  };
}

vector<Conserved> calc_fluxes
(Tessellation const& tessellation,
 vector<Primitive> const& cells,
 double dt,
 double time,
 SpatialReconstruction& interpolation,
 vector<Vector2D> const& facevelocity,
 HydroBoundaryConditions const& boundaryconditions,
 RiemannSolver const& rs,
 vector<CustomEvolution*> const& CellsEvolve,
 CustomEvolutionManager const& cem,
 vector<vector<double> > const& tracers)
{
  return serial_generate
    (FluxCalculator(interpolation,
		    tessellation,
		    cells,
		    tracers,
		    CellsEvolve,
		    cem,
		    time, dt,
		    boundaryconditions,
		    facevelocity,
		    rs));
}

void ExternalForceContribution
(const Tessellation& tess,
 const PhysicalGeometry& pg,
 const vector<ComputationalCell>& cells,
 const vector<Extensive>& fluxes,
 const vector<Vector2D>& point_velocities,
 const SourceTerm& force,
 double t,
 double dt,
 vector<Extensive>& extensives)
{
  const vector<Extensive> diff = force
    (tess,pg,cells,fluxes,point_velocities,t);
  for(size_t i=0;i<static_cast<size_t>(tess.GetPointNo());++i){
    extensives[i].mass += dt*diff[i].mass;
    extensives[i].momentum += dt*diff[i].momentum;
    extensives[i].energy += dt*diff[i].energy;
    for(std::map<std::string,double>::const_iterator it =
	  diff[i].tracers.begin();
	it!=diff[i].tracers.end(); ++it)
      extensives[i].tracers[it->first] += dt*it->second;
  }
}

vector<Vector2D> get_all_mesh_points
(Tessellation const& tess)
{
  vector<Vector2D> res(static_cast<size_t>(tess.GetPointNo()));
  for(int i=0;i<static_cast<int>(tess.GetPointNo());++i)
    res[static_cast<size_t>(i)] = tess.GetMeshPoint(i);
  return res;
}

vector<CustomEvolution*> convert_indices_to_custom_evolution
(const CustomEvolutionManager& cem,const vector<size_t>& indices)
{
  vector<CustomEvolution*> res(indices.size());
  for(size_t i=0;i<res.size();++i)
    res[i] = cem.getFunction(indices[i]);
  return res;
}

vector<Primitive> make_eos_consistent
(vector<Primitive> const& vp,
 EquationOfState const& eos)
{
  vector<Primitive> res = vp;
  for(int i=0;i<static_cast<int>(vp.size());++i)
    res[static_cast<size_t>(i)] = make_eos_consistent(vp[static_cast<size_t>(i)],eos);
  return res;
}

namespace {
  /*
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
  */

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

namespace {
  class TracerFluxCalculator: public LazyList<vector<double> >
  {
  public:

    TracerFluxCalculator(const Tessellation& tess,
			 const vector<Conserved> & fluxes,
			 const vector<CustomEvolution*>& cev,
			 const CustomEvolutionManager& cem,
			 const vector<vector<double> >& tracers,
			 const vector<Primitive>& cells,
			 const HydroBoundaryConditions& hbc,
			 const double dt,
			 const double time,
			 const SpatialReconstruction& interp,
			 const vector<Vector2D>& edge_velocities,
			 const vector<double>& lengths):
      tess_(tess), fluxes_(fluxes), cev_(cev), cem_(cem),
      tracers_(tracers), cells_(cells), hbc_(hbc), dt_(dt),
      time_(time), interp_(interp),
      edge_velocities_(edge_velocities), lengths_(lengths) {}

    size_t size(void) const
    {
      return static_cast<size_t>(tess_.GetTotalSidesNumber());
    }

    vector<double> operator[](size_t i) const
    {
      const Edge& edge = tess_.GetEdge(static_cast<int>(i));
      const double dm = fluxes_[i].Mass;
      const int n1 = edge.neighbors.second;
      const int n0 = edge.neighbors.first;
      if(hbc_.IsBoundary(edge,tess_)){
	const bool b0 = hbc_.IsGhostCell(n0,tess_);
	const bool b1 = hbc_.IsGhostCell(n1,tess_);
	assert(b0!=b1 && "One cell must be a normal cell, and the other must be a ghost cell");
	const int rci = b1 ? n0 : n1;
	return hbc_.CalcTracerFlux
	  (tess_,cells_,tracers_,dm,edge,rci,dt_,time_,interp_,
	   edge_velocities_[i]);
      }
      else{
	if(cev_[static_cast<size_t>(n0)]||cev_[static_cast<size_t>(n1)]){
	  const int ns = choose_special_cell_index(cev_,cem_,n0,n1);
	  return cev_[static_cast<size_t>(ns)]->CalcTracerFlux
	    (tess_,cells_,tracers_,dm,edge,ns,dt_,
	     time_,interp_,edge_velocities_[i]);
	}
	else
	  return apply_to_each_term
	    (interp_.interpolateTracers
	     (tess_,cells_,tracers_,dt_,edge,dm<0,InBulk,
	      edge_velocities_[i]),
	     ScalarMultiply(dt_*dm*lengths_[i]));
      }
    }

  private:
    const Tessellation& tess_;
    const vector<Conserved>& fluxes_;
    const vector<CustomEvolution*>& cev_;
    const CustomEvolutionManager& cem_;
    const vector<vector<double> >& tracers_;
    const vector<Primitive>& cells_;
    const HydroBoundaryConditions& hbc_;
    const double dt_;
    const double time_;
    const SpatialReconstruction& interp_;
    const vector<Vector2D>& edge_velocities_;
    const vector<double>& lengths_;
  };
}

void really_update_extensive_tracers
(vector<vector<double> >& extensive_tracers,
 const vector<vector<double> >& tracers,
 const vector<Primitive>& cells,
 const Tessellation& tess,
 const vector<Conserved>& fluxes,
 double time, double dt,
 const HydroBoundaryConditions& hbc,
 const SpatialReconstruction& interp,
 const vector<CustomEvolution*>& ce,
 const CustomEvolutionManager& cem,
 const vector<Vector2D>& fv,
 const vector<double>& lengths)
{
  const vector<vector<double> >& tracer_change = 
    CalcTraceChange(tracers,cells,tess,fluxes,dt,hbc,
		    interp,time,ce,cem,fv,lengths);
  UpdateTracerExtensive
    (extensive_tracers, tracer_change,
     ce, cells, tess, time);
}

vector<vector<double> > CalcTraceChange
(vector<vector<double> > const& old_trace,
 vector<Primitive> const& cells,
 Tessellation const& tess,vector<Conserved> const& fluxes,double dt,
 HydroBoundaryConditions const& hbc,
 SpatialReconstruction const& interp,
 double time,vector<CustomEvolution*> const& CellsEvolve,
 CustomEvolutionManager const& cem,
 vector<Vector2D> const& edge_velocities,
 vector<double> const& lengths)
{
  if(old_trace.empty())
    return vector<vector<double> >();

  const vector<vector<double> >& tracer_fluxes =
    serial_generate(TracerFluxCalculator
		    (tess,fluxes,CellsEvolve,cem,
		     old_trace, cells, hbc, dt, time,
		     interp, edge_velocities, lengths));
  vector<vector<double> > res(old_trace.size(),
			      vector<double>(old_trace[0].size(),0));
  for(int i=0;i<tess.GetTotalSidesNumber();++i){
    for(size_t j=0;j<res[0].size();++j){
      if(!hbc.IsGhostCell(tess.GetEdge(i).neighbors.first,tess))
	res[static_cast<size_t>(tess.GetEdge(i).neighbors.first)][j] -= tracer_fluxes[static_cast<size_t>(i)][j];
      if(!hbc.IsGhostCell(tess.GetEdge(i).neighbors.second,tess))
	res[static_cast<size_t>(tess.GetEdge(i).neighbors.second)][j] += tracer_fluxes[static_cast<size_t>(i)][j];
    }
  }
  return res;
}

vector<double> GetMaxKineticEnergy(Tessellation const& tess,vector<Primitive> const&
				   cells,vector<CustomEvolution*> const& /*customevolve*/)
{
  const int n=tess.GetPointNo();
  vector<double> res;
  res.resize(static_cast<size_t>(n));
  for(int j=0;j<n;++j)
    {
      vector<int> neightemp=tess.GetNeighbors(j);
      vector<int> neigh;
      for(size_t i=0;i<neightemp.size();++i)
	if(neightemp[i]>=0)
	  neigh.push_back(neightemp[i]);
      double e=pow(abs(cells[static_cast<size_t>(j)].Velocity-cells[static_cast<size_t>(neigh[0])].Velocity),2);
      for(int i=1;i<static_cast<int>(neigh.size());++i)
	{// This could be made much faster by writing the expression implicitly
	  e=max(e,pow(abs(cells[static_cast<size_t>(j)].Velocity-cells[static_cast<size_t>(neigh[static_cast<size_t>(i)])].Velocity),2));
	}
      res[static_cast<size_t>(j)]=0.5*e;
    }
  return res;
}

vector<double> GetForceEnergy(Tessellation const& tess,
			      vector<double> const& g)
{
  vector<double> res;
  int n=int(g.size());
  res.resize(static_cast<size_t>(n));
  for(int i=0;i<n;++i)
    res[static_cast<size_t>(i)]=g[static_cast<size_t>(i)]*tess.GetWidth(i);
  return res;
}

void FixPressure(vector<Conserved> &intensive,vector<vector<double> > const& entropy,
		 EquationOfState const& eos,vector<double> const& Ek,
		 vector<double> const& Ef,double as,double bs,vector<CustomEvolution*>
		 const& customevolve,Tessellation const& tess,//vector<Conserved> &extensive,
		 vector<char> const& shockedcells,bool densityfloor)
{
  int n=tess.GetPointNo();
  double Et,Ek2;
  double temp;
  for(int i=0;i<n;++i)
    {
      if(customevolve[static_cast<size_t>(i)]==0||customevolve[static_cast<size_t>(i)]->TimeStepRelevant())
	{
	  if (intensive[static_cast<size_t>(i)].Mass < 0)
	    continue;
	  //Make intensive
	  temp=entropy[static_cast<size_t>(i)][0]/(tess.GetVolume(i)*intensive[static_cast<size_t>(i)].Mass);
	  Ek2=0.5*pow(abs(intensive[static_cast<size_t>(i)].Momentum)/intensive[static_cast<size_t>(i)].Mass,2);
	  Et=intensive[static_cast<size_t>(i)].Energy/intensive[static_cast<size_t>(i)].Mass-Ek2;
	  if((Et<as*Ek[static_cast<size_t>(i)])||(Et<bs*Ef[static_cast<size_t>(i)]))
	    {
	      if((shockedcells[static_cast<size_t>(i)]==0)||Et<0)
		{
		  Et=eos.dp2e(intensive[static_cast<size_t>(i)].Mass,
			      eos.sd2p(temp,intensive[static_cast<size_t>(i)].Mass));
		  if(Et<0&&!densityfloor)
		    {
		      UniversalError eo("Negative thermal enegry");
		      eo.AddEntry("Cell index",i);
		      eo.AddEntry("Thermal energy",Et);
		      eo.AddEntry("ShockedStatus",shockedcells[static_cast<size_t>(i)]);
		      eo.AddEntry("Extensive entropy",entropy[static_cast<size_t>(i)][0]);
		      eo.AddEntry("The density",intensive[static_cast<size_t>(i)].Mass);
		      throw eo;
		    }
		  intensive[static_cast<size_t>(i)].Energy=intensive[static_cast<size_t>(i)].Mass*(Et+Ek2);
		  //extensive[i].Energy=tess.GetVolume(i)*intensive[i].Energy;
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
      if(neigh[static_cast<size_t>(i)]<0)
	return true;
      /*if(customevolve[neigh[i]]!=0)
	return true;*/
    }
  return false;
}

namespace {
  vector<double> scalar_mult(const vector<double>& v,
			     double s)
  {
    if(v.empty())
      return vector<double>();
    vector<double> res(v.size());
    for(size_t i=0;i<v.size();++i)
      res[i] = s*v[i];
    return res;
  }
}

namespace {
  class ExtensiveTracerCalculator: public LazyList<vector<double> >
  {
  public:

    ExtensiveTracerCalculator(const vector<vector<double> >& tracers,
			      const Tessellation& tess,
			      const vector<Primitive>& cells,
			      const PhysicalGeometry& pg):
      tracers_(tracers), tess_(tess), cells_(cells), pg_(pg) {}

    size_t size(void) const
    {
      if(tracers_.empty())
	return 0;
      else
	return static_cast<size_t>(tess_.GetPointNo());
    }

    vector<double> operator[](size_t i) const
    {
      return scalar_mult
	(tracers_[i],
	 pg_.calcVolume(serial_generate(CellEdgesGetter(tess_,static_cast<int>(i))))*
	 cells_[i].Density);
    }

  private:
    const vector<vector<double> >& tracers_;
    const Tessellation& tess_;
    const vector<Primitive>& cells_;
    const PhysicalGeometry& pg_;
  };
}

vector<vector<double> > calc_extensive_tracer
(const vector<vector<double> >& intensive_tracer,
 const Tessellation& tess,
 const vector<Primitive>& cells,
 const PhysicalGeometry& pg)
{
  return serial_generate(ExtensiveTracerCalculator(intensive_tracer,
						   tess,
						   cells,
						   pg));
}

void MakeTracerExtensive(vector<vector<double> > const &tracer,
			 Tessellation const& tess,
			 vector<Primitive> const& cells,
			 vector<vector<double> > &result)
{
  const size_t n=static_cast<size_t>(tess.GetPointNo());
  result.resize(n);
  for(size_t i=0;i<n;++i)
    result[i] = scalar_mult(tracer[i],
			    tess.GetVolume(static_cast<int>(i))*cells[i].Density);
}

namespace {
  vector<double> scalar_div(const vector<double>& v,
			    const double s)
  {
    vector<double> res(v.size());
    for(size_t i=0;i<res.size();++i)
      res[i] = v[i]/s;
    return res;
  }

  class IntensiveTracerCalculator: public LazyList<vector<double> >
  {
  public:

    IntensiveTracerCalculator(const vector<vector<double> >& extensive,const Tessellation& tess,
			      const vector<Primitive>& cells,const PhysicalGeometry& pg,
			      const vector<vector<double> > &old_trace_intensive,const vector<bool> &min_density,
			      const vector<CustomEvolution*> &cevolve) :
      extensive_(extensive),tess_(tess), cells_(cells), pg_(pg),
      old_trace_intensive_(old_trace_intensive), min_density_(min_density), cevolve_(cevolve){}

    size_t size(void) const
    {
      if(extensive_.empty())
	return 0;
      else
	return static_cast<size_t>(tess_.GetPointNo());
    }

    vector<double> operator[](size_t i) const
    {
      if (min_density_[i] && !cevolve_[i])
	return old_trace_intensive_[i];
      else
	{
	  const double mass = cells_[i].Density*
	    pg_.calcVolume(serial_generate(CellEdgesGetter(tess_, static_cast<int>(i))));
	  return scalar_div(extensive_[i], mass);
	}
    }

  private:
    const vector<vector<double> >& extensive_;
    const Tessellation& tess_;
    const vector<Primitive>& cells_;
    const PhysicalGeometry& pg_;
    const vector<vector<double> > &old_trace_intensive_;
    const vector<bool> &min_density_;
    const vector<CustomEvolution*> &cevolve_;
  };
}

void MakeTracerIntensive(vector<vector<double> > &tracer,
			 const vector<vector<double> >& extensive,
			 const Tessellation& tess,
			 const vector<Primitive>& cells,
			 const PhysicalGeometry& pg,vector<bool> const& min_density_on,vector<vector<double> > const& old_trace,
			 vector<CustomEvolution*> const& cevolve)
{
  tracer = serial_generate(IntensiveTracerCalculator(extensive,tess,cells,pg,old_trace,min_density_on,cevolve));
}

void UpdateTracerExtensive(vector<vector<double> > &tracerextensive,
			   vector<vector<double> > const& tracerchange,vector<CustomEvolution*> const&
			   CellsEvolve,vector<Primitive> const& cells,Tessellation const& tess,
			   double time)
{
  for(size_t i=0;i<tracerextensive.size();++i)
    if(CellsEvolve[i])
      tracerextensive[i]=CellsEvolve[i]->UpdateTracer
	(static_cast<int>(i),tracerextensive,tracerchange,cells,tess,time);
    else
      for(size_t j=0;j<tracerextensive[i].size();++j)
	tracerextensive[i][j]+=tracerchange[i][j];
}

void TracerResetCalc
(double alpha,SpatialDistribution const& originalD,
 SpatialDistribution const& originalP,SpatialDistribution const& originalVx,
 SpatialDistribution const& originalVy,vector<SpatialDistribution const*> const& originalTracers,vector<Primitive> &cells,
 Tessellation const& tess,vector<vector<double> > &tracer,
 int tracerindex,EquationOfState const& eos,vector<CustomEvolution*>
 const& cevolve,bool coldflows)
{
  const int n = tess.GetPointNo();
  if(n<1)
    return;
  Vector2D velocity;
  if(tracer.empty())
    return;
  if(tracerindex>=static_cast<int>(tracer[0].size())||tracerindex<0)
    throw UniversalError("Error in tracerReset, wrong dimension for tracer");
  for(int i=0;i<n;++i)
    {
      bool customforce=false;
      if(cevolve[static_cast<size_t>(i)])
	customforce=cevolve[static_cast<size_t>(i)]->ShouldForceTracerReset();
      if((tracer[static_cast<size_t>(i)][static_cast<size_t>(tracerindex)]<alpha)||customforce)
	{
	  velocity.Set(originalVx(tess.GetCellCM(i)),
		       originalVy(tess.GetCellCM(i)));
	  cells[static_cast<size_t>(i)]=CalcPrimitive(originalD(tess.GetCellCM(i)),
						      originalP(tess.GetCellCM(i)),velocity,eos);
	  if(tracer[static_cast<size_t>(i)][static_cast<size_t>(tracerindex)]<0)
	    tracer[static_cast<size_t>(i)][static_cast<size_t>(tracerindex)]=0;
	  for (size_t j = 0;j<tracer[static_cast<size_t>(i)].size();++j)
	    if (coldflows&&j == 0)
	      tracer[static_cast<size_t>(i)][j] = eos.dp2s(cells[static_cast<size_t>(i)].Density, cells[static_cast<size_t>(i)].Pressure);
	    else
	      if (static_cast<int>(j) != tracerindex)
		tracer[static_cast<size_t>(i)][j] = originalTracers[j]->operator()(tess.GetCellCM(i));
	}
    }
  return;
}

void GetPointToRemove(Tessellation const& tess,Vector2D const& point,
		      double R,vector<int> & PointToRemove,int Inner)
{
  int n=tess.GetPointNo();
  PointToRemove.clear();
  for(int i=Inner;i<n;++i)
    {
      // Check if point is completly engulfed
      bool test=true;
      vector<int> neigh=tess.GetNeighbors(i);
      for(int j=0;j<static_cast<int>(neigh.size());++j)
	if(neigh[static_cast<size_t>(j)]>=Inner)
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
    Vector2D par=edge.vertices.second-edge.vertices.first;
    if(fabs(par.x)>fabs(par.y)){
      // We are in the x direction
      if(MeshPoint.y>edge.vertices.first.y)
	MeshPoint.y-=MeshPoint.y-edge.vertices.first.y;
      else
	MeshPoint.y+=edge.vertices.first.y-MeshPoint.y;
    }
    else{
      // We are in the y direction
      if(MeshPoint.x>edge.vertices.first.x)
	MeshPoint.x-=MeshPoint.x-edge.vertices.first.x;
      else
	MeshPoint.x+=edge.vertices.first.x-MeshPoint.x;
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
  for(size_t i=0;i<edges.size();++i)
    edges[i]=tess.GetEdge(edges_loc[i]);

  // Calculate gradiant by gauss theorem
  double vx_i = cells[static_cast<size_t>(index)].Velocity.x;
  double vy_i = cells[static_cast<size_t>(index)].Velocity.y;
  Vector2D center=tess.GetMeshPoint(index);
  Vector2D c_ij,r_ij;
  for(size_t i=0;i<edges.size();i++)
    {
      if(hbc.IsBoundary(edges[i],tess))
	{
	  Primitive other2=hbc.GetBoundaryPrimitive(edges[i],tess,cells,
						    time);
	  if((edges[i].neighbors.first==-1)||(edges[i].neighbors.second==-1))
	    {
	      c_ij = CalcCentroid(edges[i])-0.5*(GetReflectedPoint
						 (tess,index,edges[i])+center);
	      r_ij = center-GetReflectedPoint(tess,index,edges[i]);
	    }
	  else
	    {
	      const int other_point = get_other_index(edges[i],index);
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
	  const int other = get_other_index(edges[i],index);
	  c_ij = CalcCentroid(edges[i])-
	    0.5*(tess.GetMeshPoint(other)+center);
	  r_ij = center-tess.GetMeshPoint(other);
	  double rij_1=1/abs(r_ij);
	  double vx_j = cells[static_cast<size_t>(other)].Velocity.x;
	  double vy_j = cells[static_cast<size_t>(other)].Velocity.y;
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

void FixAdvection(vector<Conserved>& extensive,
		  vector<Conserved> const& intensive,Tessellation const& tessold,
		  Tessellation const& tessnew,vector<Vector2D> const& facevelocity,
		  double dt,vector<Vector2D> const& /*pointvelocity*/)
{
  int n=tessold.GetTotalSidesNumber();
  int npoints=tessold.GetPointNo();
  vector<double> Rold(static_cast<size_t>(npoints)),Rnew(static_cast<size_t>(npoints));
  vector<vector<Vector2D> > pold(static_cast<size_t>(npoints)),pnew(static_cast<size_t>(npoints));
  for(int i=0;i<npoints;++i)
    {
      Rold[static_cast<size_t>(i)]=tessold.GetWidth(i);
      Rnew[static_cast<size_t>(i)]=tessnew.GetWidth(i);
      ConvexHull(pold[static_cast<size_t>(i)],&tessold,i);
      ConvexHull(pnew[static_cast<size_t>(i)],&tessnew,i);
    }

  PolygonOverlap polyoverlap;
  double eps=1e-7;
  for(int i=0;i<n;++i)
    {
      Edge const& edge=tessold.GetEdge(i);
      int n0=edge.neighbors.first;
      int n1=edge.neighbors.second;
      if(n0<0||n1<0)
	continue;
      Vector2D norm(tessold.GetMeshPoint(n1)-tessold.GetMeshPoint(n0));
      norm=norm/abs(norm);
      norm=norm*edge.GetLength();
      double dv_dt=ScalarProd(facevelocity[static_cast<size_t>(i)],norm)*dt;
      /*		vector<Vector2D> poly0,poly1;
			ConvexHull(poly0,&tessold,tessold.GetOriginalIndex(n0));
			ConvexHull(poly1,&tessnew,tessold.GetOriginalIndex(n1));
			if(n0>=npoints)
			{
			const Vector2D diff(tessold.GetMeshPoint(tessold.GetOriginalIndex(n0))
			-tessold.GetMeshPoint(n0));
			int N=static_cast<int>(poly0.size());
			for(int j=0;j<N;++j)
			poly0[j]-=diff;
			}
			if(n1>=npoints)
			{
			const Vector2D diff(tessnew.GetMeshPoint(tessnew.GetOriginalIndex(n1))
			-tessnew.GetMeshPoint(n1));
			int N=static_cast<int>(poly1.size());
			for(int j=0;j<N;++j)
			poly1[j]-=diff;
			}*/
      double real_dv1=polyoverlap.polygon_overlap_area
	    (pold[static_cast<size_t>(tessold.GetOriginalIndex(n0))],
	     pnew[static_cast<size_t>(tessold.GetOriginalIndex(n1))],
	     Rold[static_cast<size_t>(tessold.GetOriginalIndex(n0))]*eps,
	     Rnew[static_cast<size_t>(tessold.GetOriginalIndex(n1))]*eps);
      double real_dv0=polyoverlap.polygon_overlap_area
	    (pnew[static_cast<size_t>(tessold.GetOriginalIndex(n0))],
	     pold[static_cast<size_t>(tessold.GetOriginalIndex(n1))],
	     Rnew[static_cast<size_t>(tessold.GetOriginalIndex(n0))]*eps,
	     Rold[static_cast<size_t>(tessold.GetOriginalIndex(n1))]*eps);

      if(dv_dt>0)
	{
	  if(n0<npoints)
	    {
	      extensive[static_cast<size_t>(n0)]+=(real_dv0-dv_dt)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
	      extensive[static_cast<size_t>(n0)]-=real_dv1*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
	    }
	  if(n1<npoints)
	    {
	      extensive[static_cast<size_t>(n1)]+=(dv_dt-real_dv0)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
	      extensive[static_cast<size_t>(n1)]+=real_dv1*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
	    }
	}
      else
	{
	  if(n0<npoints)
	    {
	      extensive[static_cast<size_t>(n0)]-=(real_dv1+dv_dt)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
	      extensive[static_cast<size_t>(n0)]+=real_dv0*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
	    }
	  if(n1<npoints)
	    {
	      extensive[static_cast<size_t>(n1)]-=(-dv_dt-real_dv1)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
	      extensive[static_cast<size_t>(n1)]-=real_dv0*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
	    }
	}
    }
}

double determine_time_step(double hydro_time_step,
			   double external_dt,
			   double current_time,
			   double end_time)
{
  double dt = hydro_time_step;
  if(external_dt>0)
    dt = std::min(external_dt,dt);
  if(end_time>0)
    dt = std::min(end_time-current_time,dt);

#ifdef RICH_MPI
  double dt_temp = dt;
  MPI_Reduce(&dt_temp,&dt,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

  return dt;
}
