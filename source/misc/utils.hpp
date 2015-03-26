/*!  \file utils.hpp
  \brief Various useful functions
  \author Almog Yalinewich
*/

#ifndef UTILS_HPP
#define UTILS_HPP 1

#include <vector>
#include <algorithm>
#include "universal_error.hpp"
#include <cassert>

using std::vector;

/*! \brief Checks whether a number is a nan
  \param x A number
  \return True if nan false otherwise
*/
bool is_nan(double x);

/*! \brief Uniformly spaced array (like the matlab function with the same name)
  \param xl Lower bound
  \param xh Higher bound
  \param n Number of elements
  \return The uniformly spaced array
*/
vector<double> linspace(double xl, double xh, int n);

/*! \brief Uniformly spaced array (line numpy function with the same name)
  \param x_min Lower bound
  \param x_max Upper bound
  \param dx Difference between consecutive terms
  \return Uniformly spaced array
*/
vector<double> arange(double x_min, double x_max, double dx);

/*! \brief Term by term addition for vectors
  \param v1 Right argument
  \param v2 Left argument
  \return The addition of the two vectors
*/
template <class T> vector<T> operator+
(vector<T> const& v1, vector<T> const& v2)
{
  assert(v1.size()==v2.size() && "Vectors must have the same length");
  vector<T> res(v1.size());
  for(size_t i=0;i<v1.size();++i)
    res[i] = v1[i]+v2[i];
  return res;
}

/*! \brief Adds the same thing to all terms in a vector
  \param v Vector
  \param t Addition
  \return New vector with sums as terms
 */
template<class T> vector<T> operator+
(const vector<T>& v, const T& t)
{
  vector<T> res(v.size());
  for(size_t i=0, endp=v.size();i<endp;++i)
    res[i] = v[i] + t;
  return res;
}

/*! \brief Adds the same thing to all terms in a vector
  \param v Vector
  \param t Addition
  \return New vector with sums as terms
 */
template<class T> vector<T> operator+
(const T& t, const vector<T>& v)
{
  return v+t;
}

/*! \brief Adds the same thing to all terms in a vector
  \param v Vector
  \param t Addition
  \return Reference to united vector
 */
template<class T> vector<T>& operator+=
(vector<T>& v, const T& t)
{
  for(size_t i=0, endp=v.size();i<endp;++i)
    v[i] += t;
  return v;
}

/*! \brief Term by term vector subtraction
  \param v1 Left argument
  \param v2 Right argument
  \return  The subtraction between the vectors
*/
template <class T> vector<T> operator-
(vector<T> const& v1, vector<T> const& v2)
{
  if(v1.size()!=v2.size())
    throw UniversalError("Vectors must have the same length");

  vector<T> res(v1.size());
  for(size_t i=0;i<v1.size();++i)
    res[i] = v1[i]-v2[i];
  return res;
}

/*! \brief Multiplies all terms of a vector by a scalar
  \param d Scalar
  \param v vector
  \return The mulitplication of the vector
*/
template <class T> vector<T> operator*
(double d, vector<T> const& v)
{
  vector<T> res(v.size());
  for(size_t i=0;i<v.size();++i)
    res[i] = d*v[i];
  return res;
}

/*!
  \brief Linear Interpolation
  \param x The x vector, assumed sorted
  \param y y=f(x) vector
  \param xi The interpolation location
  \return f(xi)
*/
template <typename T>
T LinearInterpolation(const vector<T> &x,const vector<T> &y,T xi)
{
  typename vector<T>::const_iterator it=upper_bound(x.begin,x.end,xi);
  assert(it!=x.end && it!=x.begin &&
	 "X out of range in Linear Interpolation");
  if(*it==xi)
    return y[static_cast<size_t>(it-x.begin())];

  // Are we near the edge?
  if(it==x.end-1)
    return y[static_cast<size_t>(it-x.begin)]+(xi-*it)*
      (y[static_cast<size_t>(it-1-x.begin)]-
       y[static_cast<size_t>(it-x.begin)])/(*(it-1)-*it);
  else
    return y[static_cast<size_t>(it-x.begin)]+(xi-*it)*
      (y[static_cast<size_t>(it+1-x.begin)]-
       y[static_cast<size_t>(it-x.begin)])/(*(it+1)-*it);
}

/*! \brief Returns the minimal term in a vector
  \param v Vector
  \return The minimum of the vector
*/
double min(vector<double> const& v);

/*! \brief returns the maximal term in a vector
  \param v Vector
  \return The maximum of the vector
*/
double max(vector<double> const& v);

/*!
  \brief Removes the elements in v given by indeces
  \param v The vector to change
  \param indeces The cells to remove
*/
template <class T> void RemoveVector
(vector<T> &v,vector<int> &indeces)
{
  if(indeces.empty())
    return;
  sort(indeces.begin(),indeces.end());
  vector<T> result;
  result.reserve(v.size()-indeces.size());
  int counter=0;
  for(size_t i=0;i<static_cast<size_t>(indeces.back());++i)
    {
      if(size_t(indeces[size_t(counter)])==i)
	  ++counter;
      else
	result.push_back(v[i]);
    }
  for(size_t i=static_cast<size_t>(indeces.back())+1;i<v.size();++i)
    result.push_back(v[i]);
  v=result;
}

/*!
  \brief Returns only the values with indeces in index
  \param v The vector to check
  \param index The indeces to return
  \return The reduced vector
*/
template <class T> vector<T> VectorValues(vector<T> const&v,vector<int> const &index)
{
  if(index.empty() || v.empty())
    return vector<T> ();

  vector<T> result(index.size());
  for(size_t i=0;i<index.size();++i)
    result.at(i) = v.at(static_cast<size_t>(index.at(i)));
  return result;
}

/*!
  \brief Returns the sum of the vector
  \param v The vector to sum
  \return The sum
*/
template <class T> T VectorSum(vector<T> const&v)
{
  if(v.empty())
    return 0;
  T result = v[0];
  for(size_t i=1;i<v.size();++i)
    result += v[i];
  return result;
}

/*!
  \brief Returns a vector containing only unique elements
  \param v The input vector, must be SORTED!!
  \return The unique vector
*/
template <class T> vector<T> unique(vector<T> const& v)
{
  size_t n=v.size();
  vector<T> res;
  res.reserve(n);
  if(n==0)
    return res;
  res.push_back(v[0]);
  for (typename vector<T>::const_iterator it = v.begin() + 1; it != v.end();++it)
  if (*it == *(it - 1))
	  continue;
  else
	  res.push_back(*it);
  return res;
}

/*!
  \brief Returns a vector containing only indeces of unique elements
  \param v The input vector, must be SORTED!!
  \return The unique vector indeces
*/
template <class T> vector<int> unique_index(vector<T> const& v)
{
  if(v.empty())
    return vector<int>();

  vector<int> res;
  res.push_back(0);
  for(size_t i=1;i<v.size();++i)
    if(v[i]!=v[i-1])
      res.push_back(static_cast<int>(i));
  return res;
}

/*!
  \brief Returns only elements from vector v which are not in vector list, assumes list is sorted
  \param v The vector to change
  \param list The values to be taken out of v
  \return The new vector
*/
template <class T> vector<T> RemoveList(vector<T> const&v,vector<T> const&list)
{
  vector<T> res;
  for(size_t i=0;i<v.size();++i)
    if(!binary_search(list.begin(),list.end(),v[i]))
      res.push_back(v[i]);
  return res;
}

/*!
  \brief Removes the first occurence of val inside a vector
  \param vec The vector to change
  \param val The value to remove
*/

template <class T> void RemoveVal(vector<T> &vec,T val)
{
  for(size_t i=0;i<vec.size();++i)
    {
      if(vec[i]==val)
	{
	  vec.erase(vec.begin()+static_cast<long>(i));
	  return;
	}
    }
}

/*!
  \brief checks if val is in the vector
  \param vec The vector to check
  \param val The value to check inside the vector
  \return True if the value is in the vector, false otherwise
*/
template <class T> bool InVector(vector<T> const&vec,T val)
{
  int n=vec.size();
  for(int i=0;i<n;++i)
    if(vec[i]==val)
      return true;
  return false;
}

/*!
  \brief Returns the index of val in the vector
  \param vec The vector to check
  \param val The value to look for
  \return The index of the vector which first equals to val, throws exception if not found
*/
template <class T> int IndexInVector(vector<T> const&vec,T val)
{
  int n=vec.size();
  for(int i=0;i<n;++i)
    if(vec[i]==val)
      return i;
  throw UniversalError("Value not found in vector");
}
/*!
  \brief Rearranges the vector according to the indeces
  \param v The vector to rearrange
  \param indeces The rearrangement indeces
*/
template <class T> void ReArrangeVector(vector<T> &v,vector<int> const& indeces)
{
  const vector<T> temp=v;
  for(size_t i=0;i<v.size();++i)
    v[i]=temp[static_cast<size_t>(indeces[i])];
}
namespace
{
  template<class T> struct index_cmp
  {
    index_cmp(const T _arr) : arr(_arr) {}
    bool operator()(const size_t a, const size_t b) const
    { return arr[a] < arr[b]; }
  private:
    const T arr;
  };
}
/*! \brief Returns the indeces of a sort
  \param arr The array to sort
  \param res The indeces of the sort that is given as the output
*/
template<class T> void sort_index(const vector<T> & arr,vector<int>& res)
{
  res.resize(arr.size());
  for(size_t i=0;i<res.size();++i)
    res[i]=static_cast<int>(i);
  sort(res.begin(),res.end(),index_cmp<vector<T> >(arr));
}

/*! \brief Returns the indeces of a sort
  \param arr The array to sort
  \param res The indeces of the sort that is given as the output
*/
template<class T> void sort_index(const vector<T> & arr,vector<size_t>& res)
{
  res.resize(arr.size());
  for(size_t i=0;i<res.size();++i)
    res[i]=i;
  sort(res.begin(),res.end(),index_cmp<vector<T> >(arr));
}


/*! \brief Returns the indeces of a sort
  \param arr The array to sort
  \return The indeces of the sort that is given as the output
*/
template<class T> vector<size_t> sort_index(const vector<T> & arr)
{
  vector<size_t> res(arr.size());
  for(size_t i=0;i<res.size();++i)
    res[i]=i;
  sort(res.begin(),res.end(),index_cmp<vector<T> >(arr));
  return res;
}

/*! \brief Concatenates two vectors
  \param v1 vector
  \param v2 vector
  \return vector
*/
template<class T> vector<T> join(vector<T> const& v1,
				 vector<T> const& v2)
{
  if(v1.empty()&&v2.empty())
    return vector<T> ();
  vector<T> res(v1.size()+v2.size());
  if(!v1.empty())
    copy(v1.begin(),v1.end(),res.begin());
  if(!v2.empty())
    copy(v2.begin(),v2.end(),res.begin()+static_cast<int>(v1.size()));
  return res;
}

//! \brief BinaryOperation Binary operator
template<class T> class BinaryOperation
{
public:

  /*! \brief Evaluates the binary operation
    \param t1 First argument
    \param t2 Second argument
    \return Result
  */
  virtual T operator()(T const& t1, T const& t2) const = 0;

  virtual ~BinaryOperation(void) {}
};

/*! \brief Applies a binary operation on every pair of values from two vectors
  \param v1 First vector
  \param v2 Second vector
  \param bin_op Binary operation
  \return Vector
*/
template<class T> vector<T> binary_unite(vector<T> const& v1,
					 vector<T> const& v2,
					 BinaryOperation<T> const& bin_op)
{
  assert(v1.size()==v2.size());

  vector<T> res(v1.size());
  for(size_t i=0;i<v1.size();++i)
    res[i] = bin_op(v1[i],v2[i]);
  return res;
}

//! \brief Abstract class for an unary operation
template<class T> class UnaryOperation
{
public:

  /*! \brief Evaluate operator
    \param t Argument
    \return Result of operation
  */
  virtual T operator()(T const& t) const = 0;

  virtual ~UnaryOperation(void) {}
};

/*! \brief Applies an unary operator to all terms in a std::vector
  \param v Vector
  \param un_op Unary operator
  \return Vector
*/
template<class T> vector<T> apply_to_each_term(vector<T> const& v,
					       UnaryOperation<T> const& un_op)
{
  vector<T> res(v.size());
  for(size_t i=0;i<v.size();++i)
    res[i] = un_op(v[i]);
  return res;
}

/*! \brief Selects a member of std::pair
  \param p A pair
  \param index 0 for first member, 1 for second, error otherwise
  \return Either first or second members of pair
*/
template<class T> T pair_member(const std::pair<T,T>& p, int index)
{
  assert((0==index)||(1==index));
  return 0==index ? p.first : p.second;
}

/*! \brief Sets a member of std::pair
  \param p Pair
  \param index 0 for first, 1 for second, error otherwise
  \param val Value to be written
*/
template<class T> void set_pair_member(std::pair<T,T>& p, int index, const T& val)
{
  assert((0==index)||(1==index));
  if(0==index)
    p.first = val;
  else
    p.second = val;
}

/*! \brief Performs type casting for an entire vector
  \param source Source vector
  \return Source vector converted to new type
 */
template<class T, class S> vector<T> list_static_cast(const vector<S>& source)
{
  vector<T> res(source.size());
  for(size_t i=0;i<res.size();++i)
    res[i] = static_cast<T>(source[i]);
  return res;
}

/*! \brief Inserts all elements from one vector to the end of another
  \param subject Vector that will be modifies
  \param addendum Vector that will be added
 */
template<class T> void insert_all_to_back(vector<T>& subject, const vector<T>& addendum)
{
  if(!addendum.empty())
    subject.insert(subject.end(),addendum.begin(),addendum.end());
}

#endif // UTILS_HPP
