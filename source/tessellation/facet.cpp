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