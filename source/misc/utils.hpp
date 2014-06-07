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
	for(size_t i=0;i<(size_t)indeces[size_t(n-1)];++i)
	{
		if(size_t(indeces[size_t(counter)])==i)
		{
			++counter;
			continue;
		}
		else
			result.push_back(v[i]);
	}
	for(size_t i=(size_t)indeces[size_t(n-1)]+1;i<size_t(N);++i)
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
	if(index.empty())
		return vector<T> ();
	if(v.empty())
		return vector<T> ();
	int n=int(index.size());
	int N=int(v.size());
	vector<T> result;
	result.reserve(n);
	for(int i=0;i<n;++i)
	{
		if(index[i]<N)
			result.push_back(v[index[i]]);
		else
		{
			UniversalError eo("Vector out of bounds");
			eo.AddEntry("Tried to access index",index[i]);
			eo.AddEntry("Total length of accessed vector",N);
			throw eo;
		}
	}
	return result;
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
	return;
	//throw UniversalError("Value not found in vector");
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
	int n=(int)v.size();
	for(int i=0;i<n;++i)
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
		res[size_t(i)]=i;
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
	if(v1.empty()&&v2.empty())
		return vector<T> ();
	vector<T> res(v1.size()+v2.size());
	if(!v1.empty())
		copy(v1.begin(),v1.end(),res.begin());
	if(!v2.empty())
		copy(v2.begin(),v2.end(),res.begin()+v1.size());
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

//! \brief Ordered list whose terms are evaluated lazily
template<class T> class Index2Member
{
public:

	/*! \brief Returns the length of the list
	\return Length of the list
	*/
	virtual size_t getLength(void) const = 0;

	/*! \brief Returns a single member of the list
	\param i Index
	\return The i'th member of the list
	*/
	virtual T operator()(size_t i) const = 0;

	virtual ~Index2Member(void){}
};

/*! \brief Creates a vector from an Index2Member
\param i2m Lazily evaluated list
\return std::vector
*/
template<class T> vector<T> serial_generate(const Index2Member<T>& i2m)
{
	vector<T> res(i2m.getLength());
	for(size_t i=0, endp=res.size();i<endp;++i)
		res[i] = i2m(i);
	return res;
}

/*! \brief Multiplies all terms of std::vector with a scalar
\param v std::vector
\param s Scalar
\return Each term of v multiplied by s
*/
template<class T> vector<T> termwise_product(const vector<T>& v,
	const T& s)
{
	class Multiplier: public Index2Member<T>
	{
	public:

		Multiplier(const vector<T>& v_i,
			const T& s_i):
		v_(v_i), s_(s_i) {}

		size_t getLength(void) const
		{
			return v_.size();
		}

		T operator()(size_t i) const
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
	const size_t max_index)
{
	class Trimmer: public Index2Member<T>
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

template<class T> T pair_member(const std::pair<T,T>& p, int index)
{
  assert((0==index)||(1==index));
  return 0==index ? p.first : p.second;
}

template<class T> void set_pair_member(std::pair<T,T>& p, int index, const T& val)
{
  assert((0==index)||(1==index));
  if(0==index)
    p.first = val;
  else
    p.second = val;
}

#endif // UTILS_HPP
