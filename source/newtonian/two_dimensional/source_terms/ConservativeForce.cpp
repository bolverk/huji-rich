#include "ConservativeForce.hpp"

using std::min;
using std::sqrt;

namespace {

  Vector2D remove_parallel_component(const Vector2D& s,
				     const Vector2D& r)
  {
    return s-r*ScalarProd(r,s)/ScalarProd(r,r);
  }

  Vector2D MassFlux(const Tessellation& tess,
		    const PhysicalGeometry& pg,
		    const vector<Extensive>& fluxes,
		    int point)
  {
    Vector2D dm;
    const vector<int> edge_index=tess.GetCellEdges(point);
    const Vector2D center=tess.GetMeshPoint(point);
    for(size_t i=0;i<edge_index.size();++i)
      {
	const Edge& edge = tess.GetEdge(edge_index[i]);
	if(point==edge.neighbors.first){
	  if(edge.neighbors.second>-1){
	    dm-=pg.calcArea(edge)*fluxes[static_cast<size_t>(edge_index[i])].mass*
	      (center-tess.GetMeshPoint(edge.neighbors.second));
	  }
	  else{
	    const Vector2D p = edge.vertices.second - edge.vertices.first;
	    const Vector2D r=remove_parallel_component(center - edge.vertices.first,p);
	    dm-=2*pg.calcArea(edge)*fluxes[static_cast<size_t>(edge_index[i])].mass*r;
	  }
	}
	else
	  if(point==edge.neighbors.second){
	    if(edge.neighbors.first>-1)
	      {
		dm+=pg.calcArea(edge)*fluxes[static_cast<size_t>(edge_index[i])].mass*
		  (center-tess.GetMeshPoint(edge.neighbors.first));
	      }
	    else
	      {
		const Vector2D p = edge.vertices.second - edge.vertices.first;
		const Vector2D r = remove_parallel_component(center - edge.vertices.second, p);
		dm+=2*pg.calcArea(edge)*fluxes[static_cast<size_t>(edge_index[i])].mass*r;
	      }
	  }
	  else
	    throw UniversalError("Error in ConservativeForce MassFlux: Cell and edge are not mutual neighbors");
      }
    return dm;
  }
}

ConservativeForce::ConservativeForce(const Acceleration& acc):
  acc_(acc) {}

ConservativeForce::~ConservativeForce(void){}

namespace {
  class CellEdgesGetter: public Index2Member<Edge>
  {
  public:
    
    CellEdgesGetter(const Tessellation& tess, int n):
      tess_(tess), edge_indices_(tess.GetCellEdges(n)) {}

    size_t getLength(void) const
    {
      return edge_indices_.size();
    }

    Edge operator()(size_t i) const
    {
      return tess_.GetEdge(edge_indices_[i]);
    }

  private:
    const Tessellation& tess_;
    const vector<int> edge_indices_;
  };
}

vector<Extensive> ConservativeForce::operator()
  (const Tessellation& tess,
   const PhysicalGeometry& pg,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& fluxes,
   const vector<Vector2D>& point_velocities,
   const double t) const
{
  vector<Extensive> res(static_cast<size_t>(tess.GetPointNo()));
  for(size_t i=0;i<res.size();++i){
    const Vector2D acc = acc_
      (tess,cells,fluxes,t,static_cast<int>(i));
    const double volume = pg.calcVolume
      (serial_generate(CellEdgesGetter(tess,static_cast<int>(i))));
    res[i].mass = 0;
    res[i].momentum = volume*cells[i].density*acc;
    res[i].energy = 
      volume*cells[i].density*
      ScalarProd(point_velocities[i],acc)+
      0.5*ScalarProd(MassFlux(tess,pg,fluxes,static_cast<int>(i)),acc);
  }
  return res;
}

Acceleration::~Acceleration(void) {}
