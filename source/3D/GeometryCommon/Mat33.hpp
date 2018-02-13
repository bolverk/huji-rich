/*! \file Mat33.hpp
\brief A very simple class for a 3x3 matrix
\author Elad Steinberg
*/

#ifndef MAT3_HPP
#define MAT3_HPP 1

template <typename T>
class Mat33
{
private:
	T _data[3][3];

public:
	//! \brief Return the element at (row, col)
	inline T& at(int row, int col)
	{
		return _data[row][col];
	}

	//! \brief Return the element at (row, col)
	inline const T& at(int row, int col) const
	{
		return _data[row][col];
	}

	//! \brief Return the element at (row, col)
	inline T& operator()(int row, int col) { return at(row, col); }

	//! \brief Return the element at (row, col)
	inline const T& operator()(int row, int col) const { return at(row, col); }

	//! \brief Constructs a zeroed out matrix
	Mat33();

	//! \brief Constructs a matrix and fills all its values. An intializer_list is better, but C++ 11 isn't always supported
	Mat33(T d00, T d01, T d02,
		T d10, T d11, T d12,
		T d20, T d21, T d22);

	//! \brief Returns the matrix's determinant
	T determinant() const;

	//! \brief Return the inverse matrix
	Mat33<T> inverse()const;

	//! \brief Returns the transpose matrix
	Mat33<T> transpose()const;
};

template <typename T>
inline T Mat33<T>::determinant() const
{
	return at(0, 0)*(at(1, 1)*at(2, 2) - at(1, 2)*at(2, 1)) + at(0, 1)*(at(1, 2)*at(2, 0) - at(1, 0)*at(2, 2))
		+ at(0, 2)*(at(1, 0)*at(2, 1) - at(1, 1)*at(2, 0));
}

template<typename T>
inline Mat33<T> Mat33<T>::inverse() const
{
	Mat33<T> temp;
	double det = determinant();
	temp._data[0][0] = (_data[1][1] * _data[2][2] - _data[2][1] * _data[1][2]) / det;
	temp._data[0][1] = (_data[0][2] * _data[2][1] - _data[0][1] * _data[2][2]) / det;
	temp._data[0][2] = (_data[0][1] * _data[1][2] - _data[0][2] * _data[1][1]) / det;
	temp._data[1][0] = (_data[1][2] * _data[2][1] - _data[1][0] * _data[2][2]) / det;
	temp._data[1][1] = (_data[0][0] * _data[2][2] - _data[0][2] * _data[2][0]) / det;
	temp._data[1][2] = (_data[0][2] * _data[1][0] - _data[0][0] * _data[1][2]) / det;
	temp._data[2][0] = (_data[1][0] * _data[2][1] - _data[1][1] * _data[2][1]) / det;
	temp._data[2][1] = (_data[0][1] * _data[2][0] - _data[0][0] * _data[2][1]) / det;
	temp._data[2][2] = (_data[1][1] * _data[0][0] - _data[0][1] * _data[1][0]) / det;
	return temp;
}

template<typename T>
inline Mat33<T> Mat33<T>::transpose() const
{
	Mat33<T> res;
	for (size_t i = 0; i < 3; ++i)
		for (size_t j = 0; j < 3; ++j)
			res._data[i][j] = _data[j][i];
	return res;
}

template <typename T>
inline Mat33<T>::Mat33(T d00, T d01, T d02, 
	T d10, T d11, T d12, 
	T d20, T d21, T d22)
{
	_data[0][0] = d00;
	_data[0][1] = d01;
	_data[0][2] = d02;
	_data[1][0] = d10;
	_data[1][1] = d11;
	_data[1][2] = d12;
	_data[2][0] = d20;
	_data[2][1] = d21;
	_data[2][2] = d22;
}

template<typename T>
Mat33<T>::Mat33()
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			_data[i][j] = 0;
}

#endif // MAT33_HPP