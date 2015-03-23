#include "simple_flux_calculator.hpp"
#include "../common/hydrodynamics.hpp"
#include "../../tessellation/geometry.hpp"

SimpleFluxCalculator::SimpleFluxCalculator(const RiemannSolver& rs):
  rs_(rs) {}

namespace {
  Primitive convert_to_primitive(const ComputationalCell& cell,
				 const EquationOfState& eos)
  {
    return CalcPrimitive(cell.density,
			 cell.pressure,
			 cell.velocity,
			 eos);
  }

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

  Primitive reflect(const Primitive& p,
		    const Vector2D& axis)
  {
    return Primitive(p.Density,
		     p.Pressure,
		     Reflect(p.Velocity,axis),
		     p.Energy,
		     p.SoundSpeed);
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

  Vector2D remove_parallel_component(const Vector2D& v,
				     const Vector2D& p)
  {
    return v - p*ScalarProd(v,p)/ScalarProd(p,p);
  }

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
      res.left = reflect(res.right,res.p);
      res.n = remove_parallel_component
	(tess.GetMeshPoint(edge.neighbors.second) - edge.vertices.first,
	 res.p);
    }
    else if(!flags.second){
      res.left = convert_to_primitive
	(cells[static_cast<size_t>(edge.neighbors.first)],eos);
      res.right = reflect(res.left,res.p);
      res.n = remove_parallel_component
	(edge.vertices.first - tess.GetMeshPoint(edge.neighbors.first),
	 res.p);
    }
    else{
      const size_t left_index = static_cast<size_t>(edge.neighbors.first);
      const size_t right_index = static_cast<size_t>(edge.neighbors.second);
      res.left = convert_to_primitive(cells[left_index],eos);
      res.right = convert_to_primitive(cells[right_index],eos);
      res.velocity = ScalarProd
	(res.n,tess.CalcFaceVelocity
	 (point_velocities[left_index], 
	  point_velocities[right_index],
	  tess.GetCellCM(edge.neighbors.first), 
	  tess.GetCellCM(edge.neighbors.second),
	  calc_centroid(edge)));
      res.n = tess.GetMeshPoint(edge.neighbors.second) - 
	tess.GetMeshPoint(edge.neighbors.first);
    }
    return res;
  }
}

Conserved SimpleFluxCalculator::calcHydroFlux
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
  rpi.left.Velocity = Vector2D(Projection(rpi.left.Velocity,rpi.n),
			       Projection(rpi.left.Velocity,rpi.p));
  rpi.right.Velocity = Vector2D(Projection(rpi.right.Velocity,rpi.n),
				Projection(rpi.right.Velocity,rpi.p));
  Conserved res = rs_(rpi.left, rpi.right, rpi.velocity);
  res.Momentum = res.Momentum.x*rpi.n/abs(rpi.n)+
    res.Momentum.y*rpi.p/abs(rpi.p);
  return res;
}

vector<Extensive> SimpleFluxCalculator::operator()
(const Tessellation& tess,
 const vector<Vector2D>& point_velocities,
 const vector<ComputationalCell>& cells,
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
      res[i].tracers[it->first] = (it->second)*hydro_flux.Mass;
  }
  return res;
}
