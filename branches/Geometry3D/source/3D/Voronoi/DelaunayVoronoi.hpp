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
#include "GhostBusters.hpp"

template<typename DelaunayType, typename GhostBusterType>
class DelaunayVoronoi: public TessellationBase
{
private:
	// Make sure DelaunayType is derived from Delaunay. If you see an error on this line, it means you've instantiated DelaunayVoronoi
	// with a wrong template argument.
	BOOST_STATIC_ASSERT(boost::is_base_of<typename Delaunay, typename DelaunayType>::value);

	// Make sure GhostBusterType is derived from GhostBuster. If you see an error on this line, you've instantiated the template
	// with a wrong template argument
	BOOST_STATIC_ASSERT(boost::is_base_of<typename GhostBuster, typename GhostBusterType>::value);

public:

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

private:
	//\brief Finds a Big Tetrahedron that contains the boundray subcube, and also the 26 adjacent subcubes.
	Tetrahedron FindBigTetrahedron() const;
};

template<typename DelaunayType, typename GhostBusterType>
void DelaunayVoronoi<DelaunayType, GhostBusterType>::Update(const vector<Vector3D> &points)
{
	Tetrahedron big = FindBigTetrahedron();

	// First phase
	DelaunayType del1(points, big);
	del1.Run();

	// Find the ghost points
	GhostBusterType ghostBuster;
	set<Vector3D> ghosts = ghostBuster(del1, *_boundary);

	// Now the second phase, with the ghost points
	vector<Vector3D> allPoints(points);
	allPoints.insert(allPoints.end(), ghosts.begin(), ghosts.end());
	DelaunayType del2(allPoints, big);
	del2.Run();

}

template<typename DelaunayType, typename GhostBusterType>
Tessellation3D *DelaunayVoronoi<DelaunayType, GhostBusterType>::clone() const
{
	return new DelaunayVoronoi<DelaunayType, GhostBusterType>(*this);
}

template<typename DelaunayType, typename GhostBusterType>
vector<vector<size_t> >& DelaunayVoronoi<DelaunayType, GhostBusterType>::GetDuplicatedPoints()
{
	static vector<vector<size_t>> v;
	return v;
}

template<typename DelaunayType, typename GhostBusterType>
vector<vector<size_t> >const& DelaunayVoronoi<DelaunayType, GhostBusterType>::GetDuplicatedPoints() const
{
	static vector<vector<size_t>> v;
	return v;
}

template<typename DelaunayType, typename GhostBusterType>
size_t DelaunayVoronoi<DelaunayType, GhostBusterType>::GetTotalPointNumber() const
{
	return GetPointNo();
}

template<typename DelaunayType, typename GhostBusterType>
bool DelaunayVoronoi<DelaunayType, GhostBusterType>::IsGhostPoint(size_t index) const
{
	return false;
}
template <typename DelaunayType, typename GhostBusterType>
Tetrahedron DelaunayVoronoi<DelaunayType, GhostBusterType>::FindBigTetrahedron() const
{
	// A big tetrahedron that will contain the bounding box, as well as the 8 adjacent boundary boxes,
	// and with room to spare.

	const Vector3D &fur = _boundary->FrontUpperRight();
	const Vector3D &bll = _boundary->BackLowerLeft();
	Vector3D absFrontUpperRight(abs(fur.x), abs(fur.y), abs(fur.z));
	Vector3D absBackLowerLeft(abs(bll.x), abs(bll.y), abs(bll.z));

	absFrontUpperRight *= 1000;
	absBackLowerLeft *= -1000;

	// The top of the tetrahedron will be on the Y axis
	vector<Vector3D> tetrahedron;
	tetrahedron.push_back(Vector3D(0, absFrontUpperRight.y, 0));

	// The bottom face is parallel to the x-z plane
	double bottomY = absBackLowerLeft.y;

	// The bottom face is a triangle whose lower edge is parallel to the x axis
	double backZ = absBackLowerLeft.z;
	tetrahedron.push_back(Vector3D(absBackLowerLeft.x, bottomY, backZ));
	tetrahedron.push_back(Vector3D(absFrontUpperRight.x, bottomY, backZ));

	// The last triangle edge is on x=0
	tetrahedron.push_back(Vector3D(0, bottomY, absFrontUpperRight.z));

	return Tetrahedron(tetrahedron);
}

#endif // DELAUNAY_VORONOI_HPP