#ifndef CACHED_LAZY_LIST_HPP
#define CACHED_LAZY_LIST_HPP 1

#include <vector>
#include "lazy_list.hpp"

using std::vector;
using std::size_t;
using std::pair;

template<class T> class CachedLazyList
{
public:

  CachedLazyList(const LazyList<T>& ll):
    ll_(ll), cached_values_(ll.size(),pair<bool,T>(false,T())) {}

  size_t size(void) const
  {
    return ll_.size();
  }

  T& operator[](size_t i) const
  {
    assert(i<size());
    if(!cached_values_[i].first)
      cached_values_[i] = pair<bool,T>(true,ll_[i]);
    return cached_values_[i].second;
  }

  void reset(void) const
  {
    cached_values_.resize(ll_.size());
    for(size_t i=0;i<cached_values_.size();++i)
      cached_values_[i].first = false;
  }

  virtual ~CachedLazyList(void) {}

private:
  const LazyList<T>& ll_;
  mutable vector<pair<bool,T> > cached_values_;
};

#endif // CACHED_LAZY_LIST_HPP
