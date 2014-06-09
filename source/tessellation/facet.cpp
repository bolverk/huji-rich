#include "facet.hpp"

facet::facet():
  vertices(boost::array<int,3> ()),
  neighbors(boost::array<int,3> ()) {}


facet::facet(const facet & other):
  vertices(other.vertices),
  neighbors(other.neighbors)
{}

facet::~facet()
{
}

int facet::get_friend(int dim) const
{
	return neighbors[dim];
}

int facet::get_vertice(int dim) const
{
	return vertices[dim];
}

void facet::set_vertice(int data,int dim)
{
	vertices[dim]=data;
}

void facet::set_friend(int data,int dim)
{
	neighbors[dim]=data;
}

