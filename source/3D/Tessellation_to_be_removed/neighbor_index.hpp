#ifndef NEIGHBOR_INDEX_HPP
#define NEIGHBOR_INDEX_HPP 1

#include <vector>
#include <cassert>

using std::size_t;

class NeighborIndex
{
public:

  NeighborIndex(void);

  NeighborIndex(size_t index);

  bool isInDomain(void) const;

  size_t getIndex(void) const;

private:
  size_t index_;
  bool in_domain_;
};

#endif // NEIGHBOR_INDEX_HPP
