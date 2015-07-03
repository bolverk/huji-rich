/*! \file vector_initialiser.hpp
  \author Almog Yalinewich
  \brief Facilitates initialisation of vectors
 */

#ifndef VECTOR_INITIALISER_HPP
#define VECTOR_INITIALISER_HPP 1

#include <vector>

using std::vector;

//! \brief Class for initialising vectors
template <class T> class VectorInitialiser
{
public:

  /*! \brief Class constructor
    \param t First element
   */
  VectorInitialiser(const T& t):
    list_(1,t) {}

  VectorInitialiser(void): list_() {}

  /*! \brief Termination operator
    \return Reference to vector
   */
  const vector<T>& operator()(void) const
  {
    return list_;
  }

  /*! \brief Append operator
    \param t Next element
    \return Reference to self
   */
  VectorInitialiser& operator()(const T& t)
  {
    list_.push_back(t);
    return *this;
  }

  /*! \brief Append operator (pass by value)
    \param t Next element
    \return Reference to self
   */
  VectorInitialiser& operator[](T t)
  {
    list_.push_back(t);
    return *this;
  }

private:
  vector<T> list_;
};

#endif // VECTOR_INITIALISER_HPP
