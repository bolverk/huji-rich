//\file VectorRepository.hpp
//\brief A repository of vectors, allowing the use of vector ids instead of the entire vector
//\author Itay Zandbank

#ifndef VECTOR_REPOSITORY_HPP
#define VECTOR_REPOSITORY_HPP

#include "Vector3D.hpp"
#include <iostream>

class VectorRepository;

class VectorRef
{
private:
	size_t _id;

	VectorRef(size_t id);

public:
	const Vector3D *operator->() const;
	const Vector3D &operator*() const;

	friend class VectorRepository;
	friend std::hash<VectorRef> ;
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

//\brief Returns a VectorRef of a vector
//\remarks This operation takes O(log n), with n being the number of vectors in the repository
VectorRef GetVector(double x, double y, double z);

std::ostream& operator<<(std::ostream& output, const VectorRef &vref);

#endif