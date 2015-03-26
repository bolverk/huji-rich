#ifndef CACHED_LAZY_LIST_HPP
#define CACHED_LAZY_LIST_HPP 1

#include <vector>

using std::vector;
using std::size_t;

template<class T> class CachedLazyList
{
public:

  CachedLazyList(size_t len):
    cached_(len,false), values_(len) {}

  virtual size_t size(void) const = 0;

  virtual T calculate(const size_t i) const = 0;

  T& operator[](size_t i) const
  {
    if(cached_[i])
      return values_[i];
    else{
      cached_[i] = true;
      values_[i] = calculate(i);
      return values_[i];
    }
  }

  void reset(void) const
  {
    for(size_t i=0;i<cached_.size();++i)
      cached_[i] = false;
  }

  void resize(size_t new_size) const
  {
    cached_.resize(new_size);
    values_.resize(new_size);
    reset();
  }

  virtual ~CachedLazyList(void) {}

private:
  mutable vector<bool> cached_;
  mutable vector<T> values_;
};

#endif // CACHED_LAZY_LIST_HPP
