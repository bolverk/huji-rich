/*! \file triplet.hpp
\brief Container for three items
\author Almog Yalinewich
 */

#ifndef TRIPLET_HPP
#define TRIPLET_HPP 1

#include "universal_error.hpp"
#include <array>

//! \brief A collection of three identical references
template<class T> class TripleConstRef
{
 public:

  //! \brief Reference to first item
  const T& first;

  //! \brief Reference to second item
  const T& second;

  //! \brief Reference to third item
  const T& third;

  /*! \brief Class constructor
    \param first_i Reference to first item
    \param second_i Reference to second item
    \param third_i Reference to third item
   */
  TripleConstRef(const T& first_i,
		 const T& second_i,
		 const T& third_i):
    first(first_i),
    second(second_i),
    third(third_i) {}
};

//! \brief A collection of 3 items of the same type
template<class T> class Triplet: public std::array<T, 3>
{
 public:

  //! \brief First item
  T& first;

  //! \brief Second item
  T& second;

  //! \brief Third item
  T& third;

  /*! \brief Class constructor
    \param first_i First item
    \param second_i Second item
    \param third_i Third item
   */ 
  Triplet(const T& first_i,
	  const T& second_i,
	  const T& third_i):
    std::array<T, 3>{{first_i, second_i, third_i}},
    first((*this)[0]),
    second((*this)[1]),
    third((*this)[2]) {}

  /*! \brief Class constructor
    \param tcr References to three items
   */
  explicit Triplet(const TripleConstRef<int>& tcr):
    Triplet(tcr.first, tcr.second, tcr.third) {}

  explicit Triplet(const Triplet& other):
    Triplet(other[0], other[1], other[2]) {}

  /*! \brief Changle all items
    \param first_i First item
    \param second_i Second item
    \param third_i Third item
   */
  void set(const T& first_i,
	   const T& second_i,
	   const T& third_i)
  {
    first = first_i;
    second = second_i;
    third = third_i;
  }

  Triplet& operator=(const Triplet& other)
  {
    std::array<T,3>::operator=(other);
    return *this;
  }
};

#endif // TRIPLET_HPP
