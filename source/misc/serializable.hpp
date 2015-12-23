#ifndef SERIALIZABLE_HPP
#define SERIALIZABLE_HPP 1

#include <vector>

using std::vector;

//! \brief Base class for a serializable object
class Serializable
{
public:

  /*! \brief Returns the size of array needed to store all data
    \returns Size of array
   */
  virtual size_t getChunkSize(void) const = 0;

  /*! \brief Convert an object to an array of numbers
    \returns Array of numbers
   */
  virtual vector<double> serialize(void) const = 0;

  /*! \brief Convert an array of numbers to an object
    \param data List of numbers
   */
  virtual void unserialize(const vector<double>& data) = 0;
};

#endif // SERIALIZABLE_HPP
