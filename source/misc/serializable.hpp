#ifndef SERIALIZABLE_HPP
#define SERIALIZABLE_HPP 1

#include <vector>
#include <cassert>
#include "lazy_list.hpp"

using std::vector;
using std::size_t;

class Serializable
{
public:

  virtual size_t getChunkSize(void) const = 0;

  virtual vector<double> serialize(void) const = 0;

  virtual void unserialize(const vector<double>& data) = 0;

  virtual ~Serializable(void);
};

namespace {
  template<class T> vector<T> chunk
  (const vector<T>& source,
   size_t i_start,
   size_t i_end)
  {
    assert(i_end>i_start);
    vector<T> res(i_end-i_start);
    for(size_t i=i_start;i<i_end;++i)
      res.at(i-i_start) = source.at(i);
    return res;
  }
}

vector<double> list_serialize
(const LazyList<Serializable*>& los);

template <class S> vector<double>
list_serialize
(const LazyList<S> los)
{
  vector<double> res(los.size()*los.at(0).getChunkSize());
  size_t counter = 0;
  for(size_t i=0;i<los.size();++i){
    const vector<double> temp = los.at(i).serialize();
    BOOST_FOREACH(double x, temp){
      res.at(counter) = x;
      ++counter;
    }
  }
  return res;
}

template<class T> vector<T> list_unserialize
(const vector<double>& data,
 const T& t)
{
  assert(data.size()%t.getChunkSize()==0);
  const size_t n = data.size()/t.getChunkSize();
  vector<T> res(n);
  for(size_t i=0;i<n;++i)
    res.at(i)->unserialize
      (chunk
       (data,i*t.getChunkSize(),
	data,(i+1)*t.getChunkSize()));
  return res;
}

#endif // SERIALIZABLE_HPP
