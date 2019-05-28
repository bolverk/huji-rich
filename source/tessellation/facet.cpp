#include "facet.hpp"

facet::facet():
  vertices(0,0,0),
  neighbors(0,0,0) {}

facet::facet(const facet & other):
  vertices(other.vertices),
  neighbors(other.neighbors)
{}

facet::facet(const TripleConstRef<int>& vertices_i,
	     const TripleConstRef<int>& neighbors_i):
  vertices(vertices_i),
  neighbors(neighbors_i) {}

facet::~facet()
{
}

facet& facet::operator=(const facet& other)
{
  if(this==&other)
    return *this;
  vertices = other.vertices;
  neighbors = other.neighbors;
  return *this;
}
