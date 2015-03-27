#include "ConservativeForce.hpp"
#include "../../../misc/lazy_list.hpp"

using std::min;
using std::sqrt;

namespace {

  Vector2D remove_parallel_component(const Vector2D& s,
				     const Vector2D& r)
  {
    return s-r*ScalarProd(r,s)/ScalarProd(r,r);
  }

  Vector2D MassFlux(const Tessellation& tess,
		    const PhysicalGeometry& /*pg*/,
		    const CacheData& cd,
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
	    /*
	    dm-=pg.calcArea(edge)*fluxes[static_cast<size_t>(edge_index[i])].mass*
	      (center-tess.GetMeshPoint(edge.neighbors.second));
	    */
	    dm-=cd.areas[i]*fluxes[static_cast<size_t>(edge_index[i])].mass*
	      (center-tess.GetMeshPoint(edge.neighbors.second));
	  }
	  else{
	    const Vector2D p = edge.vertices.second - edge.vertices.first;
	    const Vector2D r=remove_parallel_component(center - edge.vertices.first,p);
	    /*
	    dm-=2*pg.calcArea(edge)*fluxes[static_cast<size_t>(edge_index[i])].mass*r;
	    */
	    dm-=2*cd.areas[i]*fluxes[static_cast<size_t>(edge_index[i])].mass*r;
	  }
	}
	else
	  if(point==edge.neighbors.second){
	    if(edge.neighbors.first>-1)
	      {
		/*
		dm+=pg.calcArea(edge)*fluxes[static_cast<size_t>(edge_index[i])].mass*
		  (center-tess.GetMeshPoint(edge.neighbors.first));
		*/
		dm+=cd.areas[static_cast<size_t>(i)]*
		  fluxes[static_cast<size_t>(edge_index[i])].mass*
		  (center-tess.GetMeshPoint(edge.neighbors.first));
	      }
	    else
	      {
		const Vector2D p = edge.vertices.second - edge.vertices.first;
		const Vector2D r = remove_parallel_component(center - edge.vertices.second, p);
		//		dm+=2*pg.calcArea(edge)*fluxes[static_cast<size_t>(edge_index[i])].mass*r;
		dm+=2*cd.areas[static_cast<size_t>(i)]*
		  fluxes[static_cast<size_t>(edge_index[i])].mass*r;
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

vector<Extensive> ConservativeForce::operator()
  (const Tessellation& tess,
   const PhysicalGeometry& pg,
   const CacheData& cd,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& fluxes,
   const vector<Vector2D>& point_velocities,
   const double t) const
{
  vector<Extensive> res(static_cast<size_t>(tess.GetPointNo()));
  for(size_t i=0;i<res.size();++i){
    if(cells[i].stickers.count("dummy")>0 && 
       cells[i].stickers.find("dummy")->second){
      res[i].mass = 0;
      res[i].energy = 0;
      res[i].momentum.x = 0;
      res[i].momentum.y = 0;
      continue;
    }
    const Vector2D acc = acc_
      (tess,cells,fluxes,t,static_cast<int>(i));
    /*
    const double volume = pg.calcVolume
      (serial_generate(CellEdgesGetter(tess,static_cast<int>(i))));
    */
    const double volume = cd.volumes[i];
    res[i].mass = 0;
    res[i].momentum = volume*cells[i].density*acc;
    res[i].energy = 
      volume*cells[i].density*
      ScalarProd(point_velocities[i],acc)+
      0.5*ScalarProd(MassFlux(tess,pg,cd,fluxes,static_cast<int>(i)),acc);
  }
  return res;
}

Acceleration::~Acceleration(void) {}
