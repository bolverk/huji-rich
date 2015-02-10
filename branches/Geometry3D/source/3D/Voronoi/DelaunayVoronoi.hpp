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
#include <boost/optional.hpp>

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
	
	//\brief Convert the DT into a VD
	//\remarks This method fills the FaceStore and Cells stuctures
	void ConvertToVoronoi(const Delaunay &del);

	typedef std::map<Vector3D, size_t> tet_map;  // A map from vertex to the tetrahedron's center (one per vertex)
	typedef std::unordered_set<size_t> tet_set;  // A set of all the unique tetrahedra (indices)

	//\brief Creates a Voronoi face, corresponding to the Delaunay Edge
	//\return The face's index in the face store, or boost::none if the face is degenerate
	//\param del The Delaunay class after its been run
	//\param vec1 The first vertex of the Delaunay edge
	//\param vec2 The second vertex of the Delaunay edge
	boost::optional<size_t> CreateFace(const Delaunay &del, const Vector3D &vec1, const Vector3D &vec2);
	
	//\brief Constructs all the cell objects
	void ConstructCells();
};

template<typename DelaunayType, typename GhostBusterType>
void DelaunayVoronoi<DelaunayType, GhostBusterType>::Update(const vector<Vector3D> &points)
{
	_meshPoints = points;
	ClearData();

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

	ConvertToVoronoi(del2);
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


template <typename DelaunayType, typename GhostBusterType>
void DelaunayVoronoi<DelaunayType, GhostBusterType>::ConvertToVoronoi(const Delaunay &del)
{
	for (size_t i = 0; i < del.NumTetrahedra(); i++)
	{
		const Tetrahedron &t = del[i];
		boost::optional<size_t> cells[4];
		GetTetrahedronIndices(t, cells);

		// Go over all the edges
		for (int i = 0; i < 3; i++)
		{
			for (int j = i + 1; j < 4; j++)
			{
				if (!cells[i].is_initialized() && !cells[j].is_initialized())  // This edge represents cells outside the boundary
					continue;

				boost::optional<size_t> faceIndex = CreateFace(del, t[i], t[j]);
				if (!faceIndex.is_initialized())
					continue; // Degenerate face

				Face &face = _faces.GetFace(faceIndex.value());

				if (cells[i].is_initialized())
					face.AddNeighbor(cells[i].value(), !cells[j].is_initialized());
				if (cells[j].is_initialized())
					face.AddNeighbor(cells[j].value(), !cells[i].is_initialized());
			}
		}
	}

	ConstructCells();
}

template <typename DelaunayType, typename GhostBusterType>
boost::optional<size_t> DelaunayVoronoi<DelaunayType, GhostBusterType>::CreateFace(const Delaunay &del, 
	const Vector3D &vec1, const Vector3D &vec2)
{
	std::set<Vector3D> vertexSet;

	vector<int> edgeNeighbors = del.EdgeNeighbors(vec1, vec2);
	for (vector<int>::iterator itNeighbor = edgeNeighbors.begin(); itNeighbor != edgeNeighbors.end(); itNeighbor++)
	{
		const Tetrahedron &neighbor = del[*itNeighbor];
		Vector3D center = neighbor.center();
		if (_boundary->inside(center))
			vertexSet.insert(center);
	}

	if (vertexSet.size() < 3) // This is a degenerate face, ignore it
		return boost::none;

	std::vector<Vector3D> vertices;
	std::copy(vertexSet.begin(), vertexSet.end(), std::back_inserter(vertices));

#ifdef _DEBUG
	if (vertices.size() > 3)  // Check that all vertices are co-planar
	{
		for (int v = 3; v < vertices.size(); v++)
		{
			Tetrahedron t(vertices[0], vertices[1], vertices[2], vertices[v]);
			BOOST_ASSERT(t.volume() == 0); // Coplanar iff the tetrahedron from the 4 vertice's volume is 0
		}
	}
#endif

	size_t faceIndex = _faces.StoreFace(vertices);
	_faces.GetFace(faceIndex).ReorderVertices();
	return faceIndex;
}

template<typename DelaunayType, typename GhostBusterType>
void DelaunayVoronoi<DelaunayType, GhostBusterType>::ConstructCells()
{
	std::vector<std::vector<size_t>> cellFaces(_meshPoints.size());

	for (size_t faceIndex = 0; faceIndex < _faces.NumFaces(); faceIndex++)
	{
		const Face& face = _faces.GetFace(faceIndex);
		BOOST_ASSERT(face.NumNeighbors() > 0); // A face with no neighbors shouldn't have been created
		
		size_t neighbor1 = face.Neighbor1().value().GetCell();
		cellFaces[neighbor1].push_back(faceIndex);

		if (face.NumNeighbors() == 2)
		{
			size_t neighbor2 = face.Neighbor2().value().GetCell();
			cellFaces[neighbor2].push_back(faceIndex);
		}
	}

	for (int i = 0; i < _meshPoints.size(); i++)
	{
		double volume = -1; // TODO: Calculate volume
		double width = -1;  // TODO: Calculate width
		Vector3D centerOfMass; // TODO: Calculate CoM
		_cells[i] = Cell(cellFaces[i], volume, width, _meshPoints[i], centerOfMass);
		_allCMs[i] = centerOfMass;
	}
}

#endif // DELAUNAY_VORONOI_HPP