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
private:
	boost::array<int,3> vertices;
	boost::array<int,3> friends;

public:

	//! \brief Defualt constructor. Sets vertices and friends to zero.
	facet();
	//! \brief Copy constructor.
	facet(const facet & other);
	//! \brief Class destructor
	~facet();
	/*! \brief Changes the index of a vertice.
	\param data The index of the new vertice. \param dim The index in the facet to change.
	*/
	void set_vertice(int data,int dim);
	//! \brief Returns a vertice. \returns The index of the vertice. \param dim The index in the facet to return.
	int get_vertice(int dim) const;
	/*! \brief Changes the index of a friend.
	\param data The index of the new friend. \param dim The index in the facet to change.
	*/
	void set_friend(int data,int dim);
	//! \brief Returns a friend. \returns The index of the friend. \param dim The index in the facet to return.
	int get_friend(int dim) const;
};

#endif //FACET_HPP