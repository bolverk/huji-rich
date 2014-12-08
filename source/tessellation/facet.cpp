#include "facet.hpp"

facet::facet():
  vertices(0,0,0),
  neighbors(0,0,0) {}

facet::facet(const facet & other):
  vertices(other.vertices),
  neighbors(other.neighbors)
{}

facet::~facet()
{
}
