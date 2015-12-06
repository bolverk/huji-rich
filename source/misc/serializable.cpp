#include "serializable.hpp"

vector<double> list_serialize
(const LazyList<Serializable*>& los)
{
  vector<double> res(los.size()*los.at(0).getChunkSize());
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
