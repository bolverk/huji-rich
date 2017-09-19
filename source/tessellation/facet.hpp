/*! \file facet.hpp
  \brief Triangulation data
  \author Elad Steinberg
 */

#ifndef FACET_HPP
#define FACET_HPP 1

#include "geometry.hpp"
#include <boost/array.hpp>
#include "../misc/triplet.hpp"

//using namespace boost;

/*!
\brief Delanauy triangle data structure. Keeps the indexes of the vertices (right hand fashion) and the neighboring facets(triangles).
\author Elad Steinberg
*/
class facet
{
public:

  //! \brief Indices of vertices
  Triplet<int> vertices;

  //! \brief Indices of neighboring facets
  Triplet<int> neighbors;

	//! \brief Defualt constructor. Sets vertices and friends to zero.
	facet();
	//! \brief Copy constructor.
	facet(const facet & other);

  /*! \brief Class constructor
    \param vertices_i Vertices indices
    \param neighbors_i Neighbor indices
   */
  facet(const TripleConstRef<int>& vertices_i,
	const TripleConstRef<int>& neighbors_i);
  
	//! \brief Class destructor
	~facet();

  facet& operator=(const facet& other);
};

#endif //FACET_HPP
