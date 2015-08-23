#include "periodic_bc.hpp"
#include "../common/hydrodynamics.hpp"
#include "../../tessellation/geometry.hpp"

PeriodicBC::PeriodicBC
(const RiemannSolver& rs):
  rs_(rs) {}

Primitive convert_to_primitive(const ComputationalCell& cell,
			       const EquationOfState& eos)
{
  /*
  return CalcPrimitive(cell.density,
		       cell.pressure,
		       cell.velocity,
		       eos);
  */
  const double energy = eos.dp2e(cell.density, cell.pressure, cell.tracers);
  const double sound_speed = eos.dp2c(cell.density, cell.pressure, cell.tracers);
  return Primitive(cell.density, cell.pressure, cell.velocity, energy, sound_speed);
}

Primitive reflect(const Primitive& p,
		  const Vector2D& axis)
{
  return Primitive(p.Density,
		   p.Pressure,
		   Reflect(p.Velocity,axis),
		   p.Energy,
		   p.SoundSpeed);
}

Vector2D remove_parallel_component(const Vector2D& v,
				   const Vector2D& p)
{
  return v - p*ScalarProd(v,p)/ScalarProd(p,p);
}

namespace {
  template<class S, class T> class Transform
  {
  public:

    virtual T operator()(const S& s) const = 0;

    virtual ~Transform(void) {}
  };

  class CellIndexValidator: public Transform<int, bool>
  {
  public:

    CellIndexValidator(const int point_no):
      point_no_(point_no) {}

    bool operator()(const int& index) const
    {
      return index>=0 && index<point_no_;
    }

  private:
    const int point_no_;
  };

  template<class S, class T> std::pair<T,T> transform_pair
  (const std::pair<S,S>& s, const Transform<S,T>& t)
  {
    return std::pair<T,T>(t(s.first),t(s.second));
  }

  class RiemannProblemInput
  {
  public:

    Primitive left;
    Primitive right;
    double velocity;
    Vector2D n;
    Vector2D p;

    RiemannProblemInput(void):
      left(), right(), velocity(0), n(), p() {}
  };

  RiemannProblemInput riemann_reduce
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const Edge& edge)
  {
    RiemannProblemInput res;
    res.p = Parallel(edge);
    const std::pair<bool,bool> flags = transform_pair
      (edge.neighbors,CellIndexValidator(tess.GetPointNo()));
    assert(flags.first || flags.second);
    if(!flags.first){
      res.right = convert_to_primitive
	(cells[static_cast<size_t>(edge.neighbors.second)],eos);
      //      res.left = reflect(res.right,res.p);
      res.left =
	convert_to_primitive
	(cells.at
	 (static_cast<size_t>
	   (tess.GetOriginalIndex
	    (edge.neighbors.first))),
	 eos);	       
      res.n = remove_parallel_component
	(tess.GetMeshPoint(edge.neighbors.second) - edge.vertices.first,
	 res.p);
    }
    else if(!flags.second){
      res.left = convert_to_primitive
	(cells.at(static_cast<size_t>(edge.neighbors.first)),eos);
      //      res.right = reflect(res.left,res.p);
      res.right = convert_to_primitive
	(cells.at
	 (static_cast<size_t>(
	  tess.GetOriginalIndex
	  (edge.neighbors.second))),
	 eos);	
      res.n = remove_parallel_component
	(edge.vertices.first - tess.GetMeshPoint(edge.neighbors.first),
	 res.p);
    }
    else{
      const size_t left_index = static_cast<size_t>(edge.neighbors.first);
      const size_t right_index = static_cast<size_t>(edge.neighbors.second);
      res.left = convert_to_primitive(cells[left_index],eos);
      res.right = convert_to_primitive(cells[right_index],eos);
      res.n = tess.GetMeshPoint(edge.neighbors.second) - 
	tess.GetMeshPoint(edge.neighbors.first);
      res.velocity = ScalarProd
	(res.n,tess.CalcFaceVelocity
	 (point_velocities[left_index], 
	  point_velocities[right_index],
	  tess.GetCellCM(edge.neighbors.first), 
	  tess.GetCellCM(edge.neighbors.second),
	  calc_centroid(edge)))/abs(res.n);
    }
    return res;
  }
}

namespace {
  Primitive rotate(const Primitive& primitive,
		   const Vector2D& n,
		   const Vector2D& p)
  {
    return Primitive(primitive.Density,
		     primitive.Pressure,
		     Vector2D(Projection(primitive.Velocity,n),
			      Projection(primitive.Velocity,p)),
		     primitive.Energy,
		     primitive.SoundSpeed);			      
  }

  Conserved rotate_back(const Conserved& c,
			const Vector2D& n,
			const Vector2D& p)
  {
    return Conserved(c.Mass,
		     c.Momentum.x*n/abs(n)+
		     c.Momentum.y*p/abs(p),
		     c.Energy);
  }
}

Conserved rotate_solve_rotate_back
(const RiemannSolver& rs,
 const Primitive& left,
 const Primitive& right,
 const double velocity,
 const Vector2D& n,
 const Vector2D& p)
{
  return rotate_back(rs(rotate(left,n,p),
			rotate(right,n,p),
			velocity),
		     n,p);
}

Conserved PeriodicBC::calcHydroFlux
(const Tessellation& tess,
 const vector<Vector2D>& point_velocities,
 const vector<ComputationalCell>& cells,
 const EquationOfState& eos,
 const size_t i) const
{
  const Edge& edge = tess.GetEdge(static_cast<int>(i));
  RiemannProblemInput rpi = riemann_reduce(tess,
					   point_velocities,
					   cells,
					   eos,
					   edge);
  return rotate_solve_rotate_back(rs_,
				  rpi.left,
				  rpi.right,
				  rpi.velocity,
				  rpi.n,
				  rpi.p);
}

namespace {
  double calc_tracer_flux(size_t i,
			  const Tessellation& tess,
			  const vector<ComputationalCell>& cells,
			  const std::string& name,
			  const Conserved& hf)
  {
    const Edge& edge = tess.GetEdge(static_cast<int>(i));
    if(hf.Mass>0 && 
       edge.neighbors.first>0 &&
       edge.neighbors.first < tess.GetPointNo())
      return hf.Mass*
	cells[static_cast<size_t>(edge.neighbors.first)].tracers.find(name)->second;
    if(hf.Mass<0 &&
       edge.neighbors.second>0 &&
       edge.neighbors.second < tess.GetPointNo())
      return hf.Mass*
	cells[static_cast<size_t>(edge.neighbors.second)].tracers.find(name)->second;
    return 0;
  }
}

vector<Extensive> PeriodicBC::operator()
(const Tessellation& tess,
 const vector<Vector2D>& point_velocities,
 const vector<ComputationalCell>& cells,
 const vector<Extensive>& /*extensives_*/,
 const CacheData& /*cd*/,
 const EquationOfState& eos,
 const double /*time*/,
 const double /*dt*/) const
{
  vector<Extensive> res(tess.getAllEdges().size());
  for(size_t i=0;i<tess.getAllEdges().size();++i){
    const Conserved hydro_flux = calcHydroFlux
      (tess, point_velocities, cells, eos, i);
    res[i].mass = hydro_flux.Mass;
    res[i].momentum = hydro_flux.Momentum;
    res[i].energy = hydro_flux.Energy;
    for(std::map<std::string,double>::const_iterator it = 
	  cells.front().tracers.begin();
	it!=cells.front().tracers.end();++it)
      res[i].tracers[it->first] =
	calc_tracer_flux(i,tess,cells,it->first,hydro_flux);
      //      res[i].tracers[it->first] = (it->second)*hydro_flux.Mass;
  }
  return res;
}
