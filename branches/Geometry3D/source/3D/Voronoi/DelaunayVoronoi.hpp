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
#include <unordered_set>
#include <unordered_map>

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

	const static double EDGE_RATIO;

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

	std::vector<VectorRef> AllPoints;

private:
	//\brief Convert the DT into a VD
	//\remarks This method fills the FaceStore and Cells stuctures
	void ConvertToVoronoi(const Delaunay &del);

	typedef std::unordered_map<VectorRef, size_t> tet_map;  // A map from vertex to the tetrahedron's center (one per vertex)
	typedef std::unordered_set<size_t> tet_set;  // A set of all the unique tetrahedra (indices)

	//\brief Creates a Voronoi face, corresponding to the Delaunay Edge
	//\return The face's index in the face store, or boost::none if the face is degenerate
	//\param del The Delaunay class after its been run
	//\param vec1 The first vertex of the Delaunay edge
	//\param vec2 The second vertex of the Delaunay edge
	boost::optional<size_t> CreateFace(const Delaunay &del, const VectorRef vec1, const VectorRef vec2);
	
	//\brief Constructs all the cell objects
	void ConstructCells();
	//\brief Split cell into tetrahedra
	//\param faces The list of cell faces
	//\return The tetrahedra, which are going to be non-overlapping and cover the entire cell's volume
	std::vector<Tetrahedron> SplitCell(const std::vector<size_t> &faces);
};

template<typename DelaunayType, typename GhostBusterType>
void DelaunayVoronoi<DelaunayType, GhostBusterType>::Update(const vector<Vector3D> &points)
{
	_meshPoints = points;
	ClearData();

	Tetrahedron big = Delaunay::CalcBigTetrahedron(*_boundary);
	vector<VectorRef> pointRefs(points.begin(), points.end());

	// First phase
	DelaunayType del1(pointRefs, big);
	del1.Run();

	// Find the ghost points
	GhostBusterType ghostBuster;
	set<VectorRef> ghosts = ghostBuster(del1, *_boundary);

	// Now the second phase, with the ghost points
	vector<VectorRef> allPointRefs(pointRefs);
	allPointRefs.insert(allPointRefs.end(), ghosts.begin(), ghosts.end());
	AllPoints = allPointRefs;
	DelaunayType del2(allPointRefs, big);
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
const double DelaunayVoronoi<DelaunayType, GhostBusterType>::EDGE_RATIO = 1e-5;

template <typename DelaunayType, typename GhostBusterType>
boost::optional<size_t> DelaunayVoronoi<DelaunayType, GhostBusterType>::CreateFace(const Delaunay &del, 
	const VectorRef vec1, const VectorRef vec2)
{
	std::vector<VectorRef> vertices;

	vector<int> edgeNeighbors = del.EdgeNeighbors(vec1, vec2);
	double firstRadius;
	for (vector<int>::iterator itNeighbor = edgeNeighbors.begin(); itNeighbor != edgeNeighbors.end(); itNeighbor++)
	{
		const Tetrahedron &neighbor = del[*itNeighbor];
		VectorRef center = neighbor.center();
		if (!_boundary->inside(*center))
			continue;

		// Check distance from the previous center
		if (vertices.size())
		{
			double dist = abs(*center - *vertices.back());
			double threshold = neighbor.radius() * EDGE_RATIO;
			if (dist < threshold)
				continue;
		}
		else
			firstRadius = neighbor.radius();  // Remember this for later, when comparing first to last

		vertices.push_back(center);
	}

	if (vertices.size() < 3) // This is a degenerate face, ignore it
		return boost::none;

	// The above loop checks the distance between neighbor i and neighbor i+1
	// we still need to check the distance between the first and last neighbots
	double dist = abs(*vertices.back() - *vertices.front());
	double threshold = firstRadius * EDGE_RATIO;
	if (dist < threshold)
		vertices.pop_back();  // Remove the last element

	if (vertices.size() < 3)  // Check for a degenerate face again
		return boost::none;

	Face face(vertices);
	size_t faceIndex = _faces.StoreFace(face.vertices);
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

	for (size_t i = 0; i < _meshPoints.size(); i++)
	{
		std::vector<Tetrahedron> tetrahedra = SplitCell(cellFaces[i]);
		double totalVolume = 0;
		Vector3D centerOfMass;

		for (size_t j = 0; j < tetrahedra.size(); j++)
			totalVolume += tetrahedra[j].volume();

		for (size_t j = 0; j < tetrahedra.size(); j++)
		{
			Vector3D weightedCenter = *tetrahedra[j].center() * tetrahedra[j].volume() / totalVolume;
			centerOfMass += weightedCenter;
		}

		_cells[i] = Cell(cellFaces[i], totalVolume, _meshPoints[i], centerOfMass);
		_allCMs[i] = centerOfMass;
	}
}

template<typename DelaunayType, typename GhostBusterType>
std::vector<Tetrahedron> DelaunayVoronoi<DelaunayType, GhostBusterType>::SplitCell(const std::vector<size_t> &faces)
{
	std::vector<Tetrahedron> tetrahedra;
	Vector3D center;
	std::unordered_set<VectorRef> considered;

	// Find the center of the cell (an average of all the vertices)
	for (std::vector<size_t>::const_iterator itFace = faces.cbegin(); itFace != faces.cend(); itFace++)
	{
		const Face& face = GetFace(*itFace);
		for (std::vector<VectorRef>::const_iterator itVertex = face.vertices.cbegin(); itVertex != face.vertices.cend(); itVertex++)
		{
			if (considered.find(*itVertex) != considered.end())  // See if we've used this vertex before
				continue;
			considered.insert(*itVertex);
			center += **itVertex;
		}
	}
	center = center / (double)considered.size();   // Average
	VectorRef centerRef(center);

	// Now create the tetrahedra, from the center to each of the faces
	for (std::vector<size_t>::const_iterator itFace = faces.cbegin(); itFace != faces.cend(); itFace++)
	{
		const Face &face = GetFace(*itFace);
		// Split the face into trianges (face[0], face[1], face[2]), (face[0], face[2], face[3]) and so on until (face[0], face[n-2], face[n-1])
		// add center to each triangle, providing the tetrahedron
		for (size_t i = 1; i < face.vertices.size() - 1; i++)
			tetrahedra.push_back(Tetrahedron(centerRef, face.vertices[0], face.vertices[i], face.vertices[i + 1]));
	}

	return tetrahedra;
}

#endif // DELAUNAY_VORONOI_HPP