#include <boost/foreach.hpp>
#include "serializable.hpp"

Serializable::~Serializable(void) {}

vector<double> list_serialize
(const vector<Serializable*>& los)
{
  Serializable* dummy = los.at(0);
  vector<double> res(los.size()*(dummy->getChunkSize()));
  size_t counter = 0;
  for(size_t i=0;i<los.size();++i){
    const vector<double> temp = los.at(i)->serialize();
    BOOST_FOREACH(double x, temp){
      res.at(counter) = x;
      ++counter;
    }
  }
  return res;
}
