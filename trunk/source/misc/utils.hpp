/*!  \file utils.hpp
  \brief Various useful functions
  \author Almog Yalinewich
*/

#ifndef UTILS_HPP
#define UTILS_HPP 1

#include <vector>
#include <algorithm>
#include "universal_error.hpp"

using namespace std;

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

/*! \brief Term by term addition for vectors
\param v1 Right argument
\param v2 Left argument
\return The addition of the two vectors
*/
template <class T> vector<T> operator+
	(vector<T> const& v1, vector<T> const& v2)
{
  if(v1.size()!=v2.size())
    throw UniversalError("Vectors must have the same length");
  vector<T> res(v1.size());
  for(int i=0;i<(int)v1.size();++i)
    res[i] = v1[i]+v2[i];
  return res;
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
	for(int i=0;i<(int)v1.size();++i)
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
	for(int i=0;i<(int)v.size();++i)
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
	typename vector<T>::const_iterator it,begin,end;
	begin=x.begin();
	end=x.end();
	it=upper_bound(begin,end,xi);
	if(it==end)
		throw UniversalError("X out of range in Linear Interpolation");
	if(it==begin)
		throw UniversalError("X out of range in Linear Interpolation");
	if(*it==xi)
		return y[(int)(it-x.begin())];
	// Are we near the edge?
	if(it==end-1)
	{
		end=it-1;
		return y[(int)(it-begin)]+(xi-*it)*(y[(int)(end-begin)]
		-y[(int)(it-begin)])/(*end-*it);
	}
	else
	{
		end=it+1;
		return y[(int)(it-begin)]+(xi-*it)*(y[(int)(end-begin)]
		-y[(int)(it-begin)])/(*end-*it);
	}
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
	int n=int(indeces.size());
	int N=int(v.size());
	vector<T> result;
	result.reserve(N-n);
	int counter=0;
	for(int i=0;i<indeces[n-1];++i)
	{
		if(indeces[counter]==i)
		{
			++counter;
			continue;
		}
		else
			result.push_back(v[i]);
	}
	for(int i=indeces[n-1]+1;i<N;++i)
		result.push_back(v[i]);
	v=result;
}


/*!
\brief Returns a vector containing only unique elements
\param v The input vector, must be SORTED!!
\return The unique vector
*/
template <class T> vector<T> unique(vector<T> const& v)
{
	int n=int(v.size());
	vector<T> res;
	if(n==0)
		return res;
	res.push_back(v[0]);
	for(int i=1;i<n;++i)
		if(v[i]==v[i-1])
			continue;
		else
			res.push_back(v[i]);
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
	int n=(int)v.size();
	vector<T> res;
	for(int i=0;i<n;++i)
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
	int n=vec.size();
	for(int i=0;i<n;++i)
	{
		if(vec[i]==val)
		{
			vec.erase(vec.begin()+i);
			return;
		}
	}
	throw UniversalError("Value not found in vector");
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
	// create copy
	vector<T> temp=v;
	// Rearrange
	for(int i=0;i<(int)v.size();++i)
		v[i]=temp[indeces[i]];
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
	for(int i=0;i<(int) res.size();++i)
		res[i]=i;
	sort(res.begin(),res.end(),index_cmp<vector<T> >(arr));
	return;
}

/*! \brief Concatenates two vectors
  \param v1 vector
  \param v2 vector
  \return vector
 */
template<class T> vector<T> join(vector<T> const& v1,
				 vector<T> const& v2)
{
  vector<T> res(v1.size()+v2.size());
  copy(v1.begin(),v1.end(),res.begin());
  copy(v2.begin(),v2.end(),res.begin()+v1.size());
  return res;
}

#endif // UTILS_HPP
