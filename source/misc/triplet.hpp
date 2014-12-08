/*! \file triplet.hpp
\brief Container for three items
\author Almog Yalinewich
 */

#ifndef TRIPLET_HPP
#define TRIPLET_HPP 1

#include "universal_error.hpp"

template<class T> class Triplet
{
 public:

  T first, second, third;

  Triplet(const T& first_i,
	  const T& second_i,
	  const T& third_i):
    first(first_i),
    second(second_i),
    third(third_i) {}

  void set(const T& first_i,
	   const T& second_i,
	   const T& third_i)
  {
    first = first_i;
    second = second_i;
    third = third_i;;
  }

  const T& operator[](size_t i) const
  {
    if(i==0)
      return first;
    else if(i==1)
      return second;
    else if(i==2)
      return third;
    else
      throw UniversalError("Access violation in triplet[]");
  }

  T& operator[](size_t i)
  {
    if(i==0)
      return first;
    else if(i==1)
      return second;
    else if(i==2)
      return third;
    else
      throw UniversalError("Access violation in triplet[]");
  }
};

template<class T> class TripleConstRef
{
 public:

  const T& first;
  const T& second;
  const T& third;

  TripleConstRef(const T& first_i,
		 const T& second_i,
		 const T& third_i):
    first(first_i),
    second(second_i),
    third(third_i) {}
};

#endif // TRIPLET_HPP
