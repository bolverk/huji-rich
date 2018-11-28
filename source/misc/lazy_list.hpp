/*! \file lazy_list.hpp
  \author Almog Yalinewich
  \brief Lazily evaluated list
 */

#ifndef LAZY_LIST_HPP
#define LAZY_LIST_HPP 1

#include <cassert>
#include <algorithm>
#include "universal_error.hpp"

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
template<class T> vector<T> serial_generate(const LazyList<T>& ll);
template<class T> vector<T> serial_generate(const LazyList<T>& ll)
{
  vector<T> res(ll.size());
  for(size_t i=0, endp=res.size();i<endp;++i)
    res[i] = ll[i];
  return res;
}

/*! \brief Multiplies all terms of std::vector with a scalar
  \param v std::vector
  \param s Scalar
  \return Each term of v multiplied by s
*/
template<class T> vector<T> termwise_product(const vector<T>& v,
					     const T& s);
template<class T> vector<T> termwise_product(const vector<T>& v,
					     const T& s)
{
  class Multiplier: public LazyList<T>
  {
  public:

    Multiplier(const vector<T>& v_i,
	       const T& s_i):
      v_(v_i), s_(s_i) {}

    size_t size(void) const
    {
      return v_.size();
    }

    T operator[](size_t i) const
    {
      return s_*v_[i];
    }

  private:
    const vector<T>& v_;
    const T& s_;
  } multiplier(v,s);

  return serial_generate(multiplier);
}

/*! \brief Trims a list and retains only a specific number of the first terms
  \param v std::vector
  \param max_index Number of terms to retain
  \return The first max_index number of terms from v
*/
template<class T> vector<T> trim_top(const vector<T>& v,
				     const size_t max_index);
template<class T> vector<T> trim_top(const vector<T>& v,
				     const size_t max_index)
{
  class Trimmer: public LazyList<T>
  {
  public:

    Trimmer(const vector<T>& v_i,
	    const size_t max_index_i):
      v_(v_i), max_index_(max_index_i) {}

    size_t getLength(void) const
    {
      return max_index_;
    }

    T operator()(size_t i) const
    {
      return v_[i];
    }

  private:
    const vector<T>& v_;
    const size_t max_index_;
  } trimmer(v,max_index);

  return serial_generate(trimmer);
}

/*!
  \brief Calculates the total number of elements in the 2d vector
  \param vec The vector to count
  \return The number of elements
*/
template<class T> int ElementNumber(vector<vector<T> > const& vec);
template<class T> int ElementNumber(vector<vector<T> > const& vec)
{
  size_t res=0;
  for(size_t i=0;i<vec.size();++i)
    res+=vec[i].size();
  return static_cast<int>(res);
}

/*!
  \brief Exchanges memebers in vec with indeces given by indeces and data given by data
  \param vec The vector to change
  \param indeces The indeces in vec to change
  \param data The data to put inside vec
*/
template<class T> void ListExchange(vector<T> &vec,vector<int> const& indeces,
				    vector<T> const& data);
template<class T> void ListExchange(vector<T> &vec,vector<int> const& indeces,
				    vector<T> const& data)
{
  if(indeces.size()!=data.size())
    {
      UniversalError eo("Matching vectors are not the same length");
      eo.AddEntry("indeces length",static_cast<double>(indeces.size()));
      eo.AddEntry("data length",static_cast<double>(data.size()));
      throw eo;
    }
  for(size_t i=0;i<indeces.size();++i)
    vec[static_cast<size_t>(indeces[i])]=data[i];
}

/*! \brief Sums terms of a lazy list
  \param i2m lazy list
  \return Sum of all terms of i2m
 */
template<class T> T lazy_sum(const LazyList<T>& i2m);
 template<class T> T lazy_sum(const LazyList<T>& i2m)
{
  T res = i2m[0];
  for(size_t i=1;i<i2m.size();++i)
    res += i2m[i];
  return res;
}

/*! \brief Finds the maximum of a lazy list
  \param i2m Lazy list
  \return Max term of lazy list
 */
 template<class T> T lazy_max(const LazyList<T>& i2m);
template<class T> T lazy_max(const LazyList<T>& i2m)
{
  T res = i2m(0);
  for(size_t i=1;i<i2m.getLength();++i)
    res = max(res,i2m(i));
  return res;
}

/*! \brief Finds the minimum of a lazy list
  \param i2m Lazy list
  \return Min term of lazy list
 */
template<class T> T lazy_min(const LazyList<T>& i2m);
 template<class T> T lazy_min(const LazyList<T>& i2m)
{
  T res = i2m[0];
  for(size_t i=1;i<i2m.size();++i)
    res = std::min(res,i2m[i]);
  return res;
}

//! \brief Converts a vector to a lazy list
template<class T> class Echo: public LazyList<T>
{
public:

  /*! \brief Class constructor
    \param v Source vector
   */
  explicit Echo(const vector<T>& v):
    v_(v) {}

  size_t size(void) const
  {
    return v_.size();
  }

  T operator[](size_t i) const
  {
    return v_[i];
  }

private:
  const vector<T>& v_;
};

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

  size_t size(void) const
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
