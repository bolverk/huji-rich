#include "facet.hpp"

facet::facet():
  vertices(boost::array<int,3> ()),
  friends(boost::array<int,3> ()) {}


facet::facet(const facet & other):
  vertices(other.vertices),
  friends(other.friends)
{}

facet::~facet()
{
}

int facet::get_friend(int dim) const
{
	return friends[dim];
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
	friends[dim]=data;
}

