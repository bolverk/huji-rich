/*! \file triplet.hpp
\brief Container for three items
\author Almog Yalinewich
 */

#ifndef TRIPLET_HPP
#define TRIPLET_HPP 1

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
