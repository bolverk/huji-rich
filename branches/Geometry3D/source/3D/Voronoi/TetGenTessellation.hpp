//\file TetGenTessellation.hpp
//\brief Uses TetGen and a ghost buster for the entire Voronoi tessallation
//\author Itay Zandbank

#ifndef TETGENTESSELLATION_HPP
#define TETGENTESSELLATION_HPP

#include "TessellationBase.hpp"
#include "GhostBusters.hpp"
#include "TetGenDelaunay.hpp"

#include <set>

template <typename GhostBusterType>
class TetGenTessellation : public TessellationBase
{
private:
	// Make sure GhostBusterType is derived from GhostBuster. If you see an error on this line, you've instantiated the template
	// with a wrong template argument
	BOOST_STATIC_ASSERT(boost::is_base_of<typename GhostBuster, typename GhostBusterType>::value);

public:
	virtual void Update(vector<Vector3D> const& points);

	virtual Tessellation3D* clone(void) const;
	virtual vector<vector<size_t> >& GetDuplicatedPoints(void);
	virtual vector<vector<size_t> >const& GetDuplicatedPoints(void)const;
	virtual size_t GetTotalPointNumber(void)const;
	virtual bool IsGhostPoint(size_t index)const;

	std::vector<VectorRef> AllPoints;

private:
	void ConvertToVoronoi(const TetGenDelaunay &del);

	//\brief Constructs all the cell objects
	void ConstructCells();
	//\brief Split cell into tetrahedra
	//\param faces The list of cell faces
	//\return The tetrahedra, which are going to be non-overlapping and cover the entire cell's volume
	std::vector<Tetrahedron> SplitCell(const std::vector<size_t> &faces);
};

template<typename GhostBusterType>
Tessellation3D *TetGenTessellation<GhostBusterType>::clone() const
{
	return new TetGenTessellation<GhostBusterType>(*this);
}

template<typename GhostBusterType>
vector<vector<size_t> >& TetGenTessellation<GhostBusterType>::GetDuplicatedPoints()
{
	static vector<vector<size_t>> v;
	return v;
}

template<typename GhostBusterType>
vector<vector<size_t> >const& TetGenTessellation<GhostBusterType>::GetDuplicatedPoints() const
{
	static vector<vector<size_t>> v;
	return v;
}

template<typename GhostBusterType>
size_t TetGenTessellation<GhostBusterType>::GetTotalPointNumber() const
{
	return GetPointNo();
}

template<typename GhostBusterType>
bool TetGenTessellation<GhostBusterType>::IsGhostPoint(size_t index) const
{
	return false;
}

template<typename GhostBusterType>
void TetGenTessellation<GhostBusterType>::Update(const vector<Vector3D> &points)
{
	_meshPoints = points;
	ClearData();

	Tetrahedron big = Delaunay::CalcBigTetrahedron(*_boundary);
	vector<VectorRef> pointRefs(points.begin(), points.end());

	// First phase
	TetGenDelaunay del1(pointRefs, big, false);
	del1.Run();

	// Find the ghost points
	GhostBusterType ghostBuster;
	set<VectorRef> ghosts = ghostBuster(del1, *_boundary);

	// Now the second phase, with the ghost points
	vector<VectorRef> allPointRefs(pointRefs);
	allPointRefs.insert(allPointRefs.end(), ghosts.begin(), ghosts.end());
	AllPoints = allPointRefs;
	TetGenDelaunay del2(allPointRefs, big, true);
	del2.Run();

	ConvertToVoronoi(del2);
}

template<typename GhostBusterType>
void TetGenTessellation<GhostBusterType>::ConvertToVoronoi(const TetGenDelaunay &del)
{
	for (size_t cellNum = 0; cellNum < _meshPoints.size(); cellNum++)
	{
		const std::vector<size_t> &tetGenFaces = del.GetVoronoiCellFaces(cellNum);
		std::vector<size_t> ourFaces;

		for (std::vector<size_t>::const_iterator itFace = tetGenFaces.begin(); itFace != tetGenFaces.end(); itFace++)
		{
			const std::vector<Vector3D> faceVertices = del.GetVoronoiFace(*itFace);
			
			size_t faceIndex = _faces.StoreFace(VectorRef::vector(faceVertices));
			ourFaces.push_back(faceIndex);
			Face face = _faces.GetFace(faceIndex);
			face.AddNeighbor(cellNum);
		}

		// Now ourFaces contains a list of faces that construct our cell - build the cell
		std::vector<Tetrahedron> tetrahedra = SplitCell(ourFaces);
		double totalVolume = 0;
		Vector3D centerOfMass;

		for (size_t j = 0; j < tetrahedra.size(); j++)
			totalVolume += tetrahedra[j].volume();

		for (size_t j = 0; j < tetrahedra.size(); j++)
		{
			Vector3D weightedCenter = *tetrahedra[j].center() * tetrahedra[j].volume() / totalVolume;
			centerOfMass += weightedCenter;
		}

		_cells[cellNum] = Cell(ourFaces, totalVolume, _meshPoints[cellNum], centerOfMass);
		_allCMs[cellNum] = centerOfMass;
	}
}

template<typename GhostBusterType>
std::vector<Tetrahedron> TetGenTessellation<GhostBusterType>::SplitCell(const std::vector<size_t> &faces)
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



#endif