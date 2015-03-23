#include "extensive.hpp"

Extensive::Extensive(void):
  mass(0),
  energy(0),
  momentum(0,0),
  tracers() {}

Extensive& Extensive::operator-=(const Extensive& diff)
{
  mass -= diff.mass;
  energy -= diff.energy;
  momentum -= diff.momentum;
  
  for(std::map<std::string,double>::iterator it=tracers.begin();
      it!=tracers.end();++it)
    it->second -= diff.tracers.find(it->first)->second;

  return *this;
}

Extensive& Extensive::operator+=(const Extensive& diff)
{
  mass += diff.mass;
  energy += diff.energy;
  momentum += diff.momentum;
  
  for(std::map<std::string,double>::iterator it=tracers.begin();
      it!=tracers.end();++it)
    it->second += diff.tracers.find(it->first)->second;

  return *this;
}

Extensive operator*(const double s,
		    const Extensive& e)
{
  Extensive res;
  res.mass = s*e.mass;
  res.energy = s*e.energy;
  res.momentum = s*e.momentum;
  res.tracers = e.tracers;
  for(std::map<std::string,double>::iterator it=res.tracers.begin();
      it!=res.tracers.end(); ++it)
    it->second *= s;
  return res;
}
