//\file DelaunayVoronoi.hpp
//\brief A Tessellation3D implementation using a Delaunay and converting to to a Voronoi Diagram
//\author Itay Zandbank
//\remarks This is a template class (the template argument being the specific Delaunay implementation),
// so everything is in this file - no cpp file for us.

#ifndef DELAUNAY_VORONOI_HPP
#define DELAUNAY_VORONOI_HPP 1

#include <boost/shared_ptr.hpp>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>
#include "Delaunay.hpp"
#include "TessellationBase.hpp"

template<typename T>
class DelaunayVoronoi: public TessellationBase
{
private:
	// Make sure T is derived from Delaunay. If you see an error on this line, it means you've instantiated DelaunayVoronoi
	// with a wrong template argument.
	BOOST_STATIC_ASSERT(boost::is_base_of<typename Delaunay, typename T>::value);

public:
	DelaunayVoronoi();

	/*!
	\brief Update the tessellation
	\param points The new positions of the mesh generating points
	*/
	virtual void Update(vector<Vector3D> const& points);


	/*!
	\brief Cloning function
	*/
	virtual Tessellation3D* clone(void) const;

	/*!
	\brief Returns the indeces of the points that where sent to other processors as ghost points (or to same cpu for single thread) ad boundary points
	\return The sent points, outer vector is the index of the outer Face and inner vector are the points sent through the face
	*/
	virtual vector<vector<size_t> >& GetDuplicatedPoints(void);
	/*!
	\brief Returns the indeces of the points that where sent to other processors as ghost points (or to same cpu for single thread) ad boundary points
	\return The sent points, outer vector is the index of the outer Face and inner vector are the points sent through the face
	*/
	virtual vector<vector<size_t> >const& GetDuplicatedPoints(void)const;
	/*!
	\brief Returns the total number of points (including ghost)
	\return The total number of points
	*/
	virtual size_t GetTotalPointNumber(void)const;

	/*!
	\brief Checks if a point is a ghost point or not
	\return True if is a ghost point, false otherwise
	*/
	virtual bool IsGhostPoint(size_t index)const;


};


#endif // DELAUNAY_VORONOI_HPP