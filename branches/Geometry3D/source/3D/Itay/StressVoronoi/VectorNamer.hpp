//\file VectorNamer.hpp
//\brief A utility class that gives string name to 3D vectors, for easier printing
//\author Itay Zandbank

#ifndef VECTORNAMER_HPP
#define VECTORNAMER_HPP

#include <unordered_map>
#include <string>
#include <sstream>
#include "GeometryCommon/VectorRepository.hpp"

class VectorNamer
{
private:
	std::unordered_map<VectorRef, std::string> _vectors;
	int _last;

public:
	VectorNamer()
	{
		_last = 0;
	}

	void SaveZero()
	{
		GetName(Vector3D(), "Z");
	}

	string GetName(VectorRef vec, string prefix = "V")
	{
		auto inMap = _vectors.find(vec);
		if (inMap != _vectors.end())
			return inMap->second;

		stringstream strm;
		strm << prefix << _last++;
		_vectors[vec] = strm.str();
		return strm.str();
	}

	std::unordered_map<VectorRef, string>::const_iterator begin() const { return _vectors.cbegin(); }
	std::unordered_map<VectorRef, string>::const_iterator end() const { return _vectors.cend(); }
};

#endif