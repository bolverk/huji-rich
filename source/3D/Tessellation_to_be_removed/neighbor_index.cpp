#include "neighbor_index.hpp"

NeighborIndex::NeighborIndex(void):
  index_(0), in_domain_(false) {}

NeighborIndex::NeighborIndex(size_t index):
  index_(index), in_domain_(true) {}

bool NeighborIndex::isInDomain(void) const
{
  return in_domain_;
}

size_t NeighborIndex::getIndex(void) const
{
  assert(in_domain_);
  return index_;
}
