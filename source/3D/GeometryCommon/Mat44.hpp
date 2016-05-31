/*! \file Mat44.hpp
\brief A very simple class for a 4x4 matrix
\author Itay Zandbank
*/

#ifndef MAT44_HPP
#define MAT44_HPP 1

template <typename T>
class Mat44
{
private:
	T _data[4][4];

public:
	// \brief Return the element at (row, col)
	inline T& at(int row, int col)
	{
		return _data[row][col];
	}

	// \brief Return the element at (row, col)
	inline const T& at(int row, int col) const
	{
		return _data[row][col];
	}

	// \brief Return the element at (row, col)
	inline T& operator()(int row, int col) { return at(row, col); }

	// \brief Return the element at (row, col)
	inline const T& operator()(int row, int col) const { return at(row, col); }

	// \brief Constructs a zeroed out matrix
	Mat44();

	// \brief Constructs a matrix and fills all its values. An intializer_list is better, but C++ 11 isn't always supported
	Mat44(T d00, T d01, T d02, T d03,
		T d10, T d11, T d12, T d13,
		T d20, T d21, T d22, T d23,
		T d30, T d31, T d32, T d33);

	// \brief Returns the matrix's determinant
	T determinant() const;
};

template <typename T>
inline T Mat44<T>::determinant() const
{
	// Code adapted taken from here: http://stackoverflow.com/a/2937973/871910
	return
		at(0, 3) * at(1, 2) * at(2, 1) * at(3, 0) - at(0, 2) * at(1, 3) * at(2, 1) * at(3, 0) -
		at(0, 3) * at(1, 1) * at(2, 2) * at(3, 0) + at(0, 1) * at(1, 3) * at(2, 2) * at(3, 0) +
		at(0, 2) * at(1, 1) * at(2, 3) * at(3, 0) - at(0, 1) * at(1, 2) * at(2, 3) * at(3, 0) -
		at(0, 3) * at(1, 2) * at(2, 0) * at(3, 1) + at(0, 2) * at(1, 3) * at(2, 0) * at(3, 1) +
		at(0, 3) * at(1, 0) * at(2, 2) * at(3, 1) - at(0, 0) * at(1, 3) * at(2, 2) * at(3, 1) -
		at(0, 2) * at(1, 0) * at(2, 3) * at(3, 1) + at(0, 0) * at(1, 2) * at(2, 3) * at(3, 1) +
		at(0, 3) * at(1, 1) * at(2, 0) * at(3, 2) - at(0, 1) * at(1, 3) * at(2, 0) * at(3, 2) -
		at(0, 3) * at(1, 0) * at(2, 1) * at(3, 2) + at(0, 0) * at(1, 3) * at(2, 1) * at(3, 2) +
		at(0, 1) * at(1, 0) * at(2, 3) * at(3, 2) - at(0, 0) * at(1, 1) * at(2, 3) * at(3, 2) -
		at(0, 2) * at(1, 1) * at(2, 0) * at(3, 3) + at(0, 1) * at(1, 2) * at(2, 0) * at(3, 3) +
		at(0, 2) * at(1, 0) * at(2, 1) * at(3, 3) - at(0, 0) * at(1, 2) * at(2, 1) * at(3, 3) -
		at(0, 1) * at(1, 0) * at(2, 2) * at(3, 3) + at(0, 0) * at(1, 1) * at(2, 2) * at(3, 3);
}

template <typename T>
inline Mat44<T>::Mat44(T d00, T d01, T d02, T d03,
	T d10, T d11, T d12, T d13,
	T d20, T d21, T d22, T d23,
	T d30, T d31, T d32, T d33)
{
	_data[0][0] = d00;
	_data[0][1] = d01;
	_data[0][2] = d02;
	_data[0][3] = d03;
	_data[1][0] = d10;
	_data[1][1] = d11;
	_data[1][2] = d12;
	_data[1][3] = d13;
	_data[2][0] = d20;
	_data[2][1] = d21;
	_data[2][2] = d22;
	_data[2][3] = d23;
	_data[3][0] = d30;
	_data[3][1] = d31;
	_data[3][2] = d32;
	_data[3][3] = d33;
}

template<typename T>
Mat44<T>::Mat44()
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			_data[i][j] = 0;
}

#endif // MAT44_HPP