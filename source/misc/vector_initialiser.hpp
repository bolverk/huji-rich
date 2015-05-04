/*! \file vector_initialiser.hpp
  \author Almog Yalinewich
  \brief Facilitates initialisation of vectors
 */

#ifndef VECTOR_INITIALISER_HPP
#define VECTOR_INITIALISER_HPP 1

template <class T> class VectorInitialiser
{
public:

  VectorInitialiser(const T& t):
    list_(1,t) {}

  const vector<T>& operator()(void) const
  {
    return list_;
  }

  VectorInitialiser& operator()(const T& t)
  {
    list_.push_back(t);
    return *this;
  }

private:
  vector<T> list_;
};

#endif // VECTOR_INITIALISER_HPP
