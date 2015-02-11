//\file VectorRepository.cpp
//\brief Implementation of the Vector Repository
//\author Itay Zandbank

#include "Vector3D.hpp"
#include "VectorRepository.hpp"
#include "../Utilities/assert.hpp"

#include <map>

using namespace std;

class VectorRepository
{
private:
	typedef map<Vector3D, size_t> map_type;
	map_type _vectorMap;
	vector<Vector3D> _repository;

public:
	VectorRef GetVectorRef(const Vector3D &vec);
	const Vector3D &GetVector(size_t index) const;
};

VectorRef VectorRepository::GetVectorRef(const Vector3D &vec)
{
	map_type::iterator it = _vectorMap.find(vec);
	if (it != _vectorMap.end())
		return VectorRef(it->second);

	_repository.push_back(vec);
	size_t index = _repository.size() - 1;
	_vectorMap[vec] = index;
	return VectorRef(index);
}

const Vector3D& VectorRepository::GetVector(size_t index) const
{
	BOOST_ASSERT(index >= 0 && index < _repository.size());
	return _repository[index];
}

static VectorRepository theRepository;

VectorRef::VectorRef(size_t id) : _id(id) { }

const Vector3D &VectorRef::operator*() const
{
	return theRepository.GetVector(_id);
}

const Vector3D *VectorRef::operator->() const
{
	return &theRepository.GetVector(_id);
}

ostream& operator<<(ostream& output, const VectorRef &vref)
{
	output << *vref;
	return output;
}
