#include "ConservativeForce.hpp"
#include "../../../misc/lazy_list.hpp"

using std::min;
using std::sqrt;

ConservativeForce::ConservativeForce(const Acceleration& acc):
  acc_(acc) {}

ConservativeForce::~ConservativeForce(void){}

vector<Extensive> ConservativeForce::operator()
  (const Tessellation& tess,
   const PhysicalGeometry& /*pg*/,
   const CacheData& cd,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& fluxes,
   const vector<Vector2D>& /*point_velocities*/,
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
    res[i].energy = volume*cells[i].density*ScalarProd(acc,cells[i].velocity);
    /*
    const Vector2D mass_flux = MassFlux(tess,pg,cd,fluxes,static_cast<int>(i));
    res[i].energy = 
      volume*cells[i].density*
      ScalarProd(point_velocities[i],acc)+
      0.5*ScalarProd(mass_flux,acc);
    */
      //      0.5*ScalarProd(MassFlux(tess,pg,cd,fluxes,static_cast<int>(i)),acc);
  }
  return res;
}

Acceleration::~Acceleration(void) {}
