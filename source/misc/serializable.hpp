#ifndef SERIALIZABLE_HPP
#define SERIALIZABLE_HPP 1

#include <vector>

using std::vector;

class Serializable
{
public:

  virtual size_t getChunkSize(void) const = 0;

  virtual vector<double> serialize(void) const = 0;

  virtual void unserialize(const vector<double>& data) = 0;
};

#endif // SERIALIZABLE_HPP
