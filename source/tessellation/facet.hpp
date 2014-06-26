/*! \file facet.hpp
  \brief Triangulation data
  \author Elad Steinberg
 */

#ifndef FACET_HPP
#define FACET_HPP 1

#include "geometry.hpp"
#include <boost/array.hpp>

using namespace boost;
/*! 
\brief Delanauy triangle data structure. Keeps the indexes of the vertices (right hand fashion) and the neighboring facets(triangles).
\author Elad Steinberg
*/
class facet
{
public:

  //! \brief Indices of vertices
  boost::array<int,3> vertices;

  //! \brief Indices of neighboring facets
  boost::array<int,3> neighbors;

	//! \brief Defualt constructor. Sets vertices and friends to zero.
	facet();
	//! \brief Copy constructor.
	facet(const facet & other);
	//! \brief Class destructor
	~facet();
};

#endif //FACET_HPP
