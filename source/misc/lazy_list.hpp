/*! \file lazy_list.hpp
  \author Almog Yalinewich
  \brief Lazily evaluated list
 */

#ifndef LAZY_LIST_HPP
#define LAZY_LIST_HPP 1

#include <cassert>
#include "universal_error.hpp"
#include <algorithm>

using std::vector;
using std::size_t;

//! \brief Ordered list whose terms are evaluated lazily
template<class T> class LazyList
{
public:

  /*! \brief Returns the length of the list
    \return Length of the list
  */
  virtual size_t size(void) const = 0;

  /*! \brief Returns a single member of the list
    \param i Index
    \return The i'th member of the list
  */
  virtual T operator[](const size_t i) const = 0;

  /*! \brief bound checked access function
    \param i Index
    \return i'th term
   */
  T at(const size_t i) const
  {
    assert(i<size());
    return (*this)[i];
  }

  virtual ~LazyList(void) {}
};

/*! \brief Creates a vector from an LazyList
  \param ll Lazily evaluated list
  \return std::vector
*/
template<class T> vector<T> serial_generate(const LazyList<T>& ll)
{
  vector<T> res(ll.size());
  for(size_t i=0, endp=res.size();i<endp;++i)
    res[i] = ll[i];
  return res;
}

//! \brief Creates a contiguous chunk of a lazy list
template<class T> class ContiguousChunk: public LazyList<T>
{
public:

  /*! \brief Class constructor
    \param i2m Source list
    \param low Lower boundary of chunk
    \param high Higher boundary of chunk
   */
  ContiguousChunk(const LazyList<T>& i2m,
		  size_t low,
		  size_t high):
    i2m_(i2m), low_(low), high_(high) 
  {
    assert(high_>low_);
  }

  size_t size(void) const override
  {
    return high_-low_;
  }

  T operator[](size_t i) const
  {
    return i2m_(i+low_);
  }

private:
  const LazyList<T>& i2m_;
  const size_t low_;
  const size_t high_;
};

#endif // LAZY_LIST_HPP
