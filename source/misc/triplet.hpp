/*! \file triplet.hpp
\brief Container for three items
\author Almog Yalinewich
 */

#ifndef TRIPLET_HPP
#define TRIPLET_HPP 1

#include "universal_error.hpp"

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
template<class T> class Triplet
{
 public:

  //! \brief First item
  T first;

  //! \brief Second item
  T second;

  //! \brief Third item
  T third;

  /*! \brief Class constructor
    \param first_i First item
    \param second_i Second item
    \param third_i Third item
   */ 
  Triplet(const T& first_i,
	  const T& second_i,
	  const T& third_i):
    first(first_i),
    second(second_i),
    third(third_i) {}

  /*! \brief Class constructor
    \param tcr References to three items
   */
  explicit Triplet(const TripleConstRef<int>& tcr):
    first(tcr.first),
    second(tcr.second),
    third(tcr.third) {}

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
    third = third_i;;
  }

  /*! \brief Random access operator
    \param i index
    \return Reference to appropriate item
   */
  const T& operator[](size_t i) const
  {
    assert(i<3);
    static const T Triplet<T>::* temp [3] = {&Triplet<T>::first,
					    &Triplet<T>::second,
					    &Triplet<T>::third};
    return this->*temp[i];
  }

  /*! \brief Random access operator
    \param i index
    \return Reference to appropriate item
   */
  T& operator[](size_t i)
  {
    assert(i<3);
    static T Triplet<T>::* temp [] = {&Triplet<T>::first,
					    &Triplet<T>::second,
					    &Triplet<T>::third};
    return this->*temp[i];
  }
};

#endif // TRIPLET_HPP
