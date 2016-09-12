#ifndef UTILS_HPP
#define UTILS_HPP 1

namespace
{
	template<class T> struct index_cmp
	{
		explicit index_cmp(const T _arr) : arr(_arr) {}
		bool operator()(const size_t a, const size_t b) const
		{
			return arr[a] < arr[b];
		}
	private:
		const T arr;
	};
}


/*! \brief Returns the indeces of a sort
\param arr The array to sort
\return The indeces of the sort that is given as the output
*/
template<class T> vector<size_t> sort_index(const vector<T> & arr)
{
	vector<size_t> res(arr.size());
	for (size_t i = 0; i < res.size(); ++i)
		res[i] = i;
	sort(res.begin(), res.end(), index_cmp<vector<T> >(arr));
	return res;
}


/*!
\brief Returns a vector containing only unique elements
\param v The input vector, must be SORTED!!
\return The unique vector
*/
template <class T> vector<T> unique(vector<T> const& v)
{
	size_t n = v.size();
	vector<T> res;
	res.reserve(n);
	if (n == 0)
		return res;
	res.push_back(v[0]);
	for (typename vector<T>::const_iterator it = v.begin() + 1; it != v.end(); ++it)
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
	if (v.empty())
		return vector<int>();

	vector<int> res;
	res.push_back(0);
	for (size_t i = 1; i < v.size(); ++i)
		if (v[i] != v[i - 1])
			res.push_back(static_cast<int>(i));
	return res;
}

/*!
\brief Rearranges the vector according to the indeces
\param v The vector to rearrange
\param indeces The rearrangement indeces
*/
template <class T> void ReArrangeVector(vector<T> &v, vector<size_t> const& indeces)
{
	const vector<T> temp = v;
	for (size_t i = 0; i < v.size(); ++i)
		v[i] = temp[indeces[i]];
}

/*! \brief Checks whether a number is a nan
\param x A number
\return True if nan false otherwise
*/
bool is_nan(double x)
{
	int b1 = (x >= 0);
	int b2 = (x<0);
	return (1 != b1 + b2);
}


#endif //UTILS_HPP