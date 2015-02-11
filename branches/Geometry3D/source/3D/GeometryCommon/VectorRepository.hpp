//\file VectorRepository.hpp
//\brief A repository of vectors, allowing the use of vector ids instead of the entire vector
//\author Itay Zandbank

#ifndef VECTOR_REPOSITORY_HPP
#define VECTOR_REPOSITORY_HPP

#include "Vector3D.hpp"
#include <iostream>

class VectorRef
{
private:
	size_t _id;

public:
	VectorRef(const Vector3D &vector);
	VectorRef();

	const Vector3D *operator->() const;
	const Vector3D &operator*() const;

	friend std::hash<VectorRef> ;
	friend bool operator==(const VectorRef &v1, const VectorRef &v2);
};

namespace std
{
	template<>
	struct hash<VectorRef>
	{
		typedef VectorRef argument_type;
		typedef size_t result_type;

		result_type operator()(const argument_type &ref) const
		{
			return ref._id;
		}
	};
}

std::ostream& operator<<(std::ostream& output, const VectorRef &vref);
bool operator==(const VectorRef &v1, const VectorRef &v2);
bool operator<(const VectorRef &v1, const VectorRef &v2);

#endif