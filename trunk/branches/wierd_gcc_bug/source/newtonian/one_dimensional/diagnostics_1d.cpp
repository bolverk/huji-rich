#include <fstream>
#include "diagnostics_1d.hpp"
#include "../../misc/simple_io.hpp"
#include "../../misc/universal_error.hpp"

double cell_property(hdsim1D const& sim,
		     int i,
		     string const& property)
{
  if("center"==property)
    return sim.GetCellCenter(i);
  else if("density"==property)
    return sim.GetCell(i).Density;
  else if("pressure"==property)
    return sim.GetCell(i).Pressure;
  else if("xvelocity"==property)
    return sim.GetCell(i).Velocity.x;
  else if("yvelocity"==property)
    return sim.GetCell(i).Velocity.y;
  else
    throw UniversalError("Unknown property "+property);
}

vector<double> cells_property(hdsim1D const& sim,
			      string const& property)
{
  vector<double> res(sim.GetCellNo(),0);
  for(int i=0;i<sim.GetCellNo();++i){
    res[i] = cell_property(sim,i,property);
  }
  return res;
}

void write_cells_property(hdsim1D const& sim,
			  string const& property,
			  string const& fname,
			  int prec)
{
  vector<double> temp = cells_property(sim,property);
  write_vector(temp,fname,prec);
}

Timer::Timer(double t_next, double dt):
  t_next_(t_next), dt_(dt), cycle_(0) {}

bool Timer::isTime(double t)
{
  if(t<t_next_)
    return false;
  else{
    t_next_ += dt_;
    ++cycle_;
    return true;
  }
}

int Timer::getCycle(void) const
{
  return cycle_;
}
