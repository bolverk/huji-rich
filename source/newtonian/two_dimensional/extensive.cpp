#include "extensive.hpp"
#include "../../misc/utils.hpp"

Extensive::Extensive(void):
  mass(0),
  energy(0),
  momentum(0,0),
  tracers() {}

Extensive::Extensive(boost::container::flat_map<std::string, double> const& Tracers):mass(0),
energy(0),
momentum(0, 0),
tracers() 
{
	for (boost::container::flat_map<std::string, double>::const_iterator it = Tracers.begin();
		it != Tracers.end(); ++it)
		tracers[it->first] = 0;
}

Extensive& Extensive::operator-=(const Extensive& diff)
{
  mass -= diff.mass;
  energy -= diff.energy;
  momentum -= diff.momentum;
  assert(diff.tracers.size() == this->tracers.size());
  boost::container::flat_map<std::string, double>::iterator it2 = this->tracers.begin();
  for (boost::container::flat_map<std::string, double>::const_iterator it = diff.tracers.begin();
  it != diff.tracers.end(); ++it, ++it2)
	  it2->second -= it->second;
  return *this;
}

Extensive& Extensive::operator=(const Extensive& origin)
{
  mass = origin.mass;
  energy = origin.energy;
  momentum = origin.momentum;
  tracers = origin.tracers;
  return *this;
}

void ReplaceExtensive(Extensive &toreplace, Extensive const& other)
{
	toreplace.mass = other.mass;
	toreplace.energy = other.energy;
	toreplace.momentum = other.momentum;
	assert(other.tracers.size() == toreplace.tracers.size());
	boost::container::flat_map<std::string, double>::iterator it2 = toreplace.tracers.begin();
	for (boost::container::flat_map<std::string, double>::const_iterator it = other.tracers.begin();
	it != other.tracers.end(); ++it, ++it2)
		it2->second = it->second;
}

Extensive& Extensive::operator+=(const Extensive& diff)
{
  mass += diff.mass;
  energy += diff.energy;
  momentum += diff.momentum;
  assert(diff.tracers.size() == this->tracers.size());
  boost::container::flat_map<std::string, double>::iterator it2 = this->tracers.begin();
  for (boost::container::flat_map<std::string, double>::const_iterator it = diff.tracers.begin();
  it != diff.tracers.end(); ++it, ++it2)
	  it2->second += it->second;
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
  for(boost::container::flat_map<std::string,double>::iterator it=res.tracers.begin();
      it!=res.tracers.end(); ++it)
    it->second *= s;
  return res;
}

Extensive operator+(const Extensive& e1,
		    const Extensive& e2)
{
  Extensive res;
  res.mass = e1.mass + e2.mass;
  res.energy = e1.energy + e2.energy;
  res.momentum = e1.momentum + e2.momentum;
  for(boost::container::flat_map<std::string,double>::const_iterator it=e1.tracers.begin();
      it!=e1.tracers.end();++it)
    res.tracers[it->first] = it->second + safe_retrieve(e2.tracers,it->first);
  return res;
}

Extensive& Extensive::operator*=(const double scalar)
{
	mass *=scalar;
	energy *=scalar;
	momentum *=scalar;
	for (size_t i = 0; i < tracers.size(); ++i)
		(tracers.begin() + static_cast<int>(i))->second *= scalar;
	return *this;
}

Extensive operator-(const Extensive& e1,
		    const Extensive& e2)
{
  return e1+(-1)*e2;
}
		   
#ifdef RICH_MPI
size_t Extensive::getChunkSize(void) const
{
  return 4 + tracers.size();
}

vector<double> Extensive::serialize(void) const
{
  vector<double> res(getChunkSize());
  res.at(0) = mass;
  res.at(1) = energy;
  res.at(2) = momentum.x;
  res.at(3) = momentum.y;
  size_t counter = 4;
  for(boost::container::flat_map<string,double>::const_iterator
	it=tracers.begin();
      it!=tracers.end();
      ++it){
    res.at(counter) = it->second;
    ++counter;
  }
  assert(counter==res.size());
  return res;
}

void Extensive::unserialize
(const vector<double>& data)
{
  assert(data.size()==getChunkSize());
  mass = data.at(0);
  energy = data.at(1);
  momentum.x = data.at(2);
  momentum.y = data.at(3);
  size_t counter = 4;
  for(boost::container::flat_map<string,double>::iterator
	it=tracers.begin();
      it!=tracers.end();
      ++it){
    it->second = data.at(counter);
    ++counter;
  }
  assert(data.size()==counter);
}
#endif // RICH_MPI
