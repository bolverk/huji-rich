//\file TetGenTessellation.hpp
//\brief Uses TetGen and a ghost buster for the entire Voronoi tessallation
//\author Itay Zandbank

#ifndef TETGENTESSELLATION_HPP
#define TETGENTESSELLATION_HPP

#include "TessellationBase.hpp"
#include "GhostBusters.hpp"
#include "TetGenDelaunay.hpp"
#include "../GeometryCommon/CellCalculations.hpp"

#include <unordered_set>

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
	const static double EDGE_RATIO;

private:
	void ConvertToVoronoi(const TetGenDelaunay &del);

	// TetGenDelaunay creates Voronoi cells for all the ghost points, too. We just need the cells of the original
	// points, which are the first _meshPoints.size() cells, and the faces that  are part of these cells.
	std::vector<boost::optional<Face>> _allFaces;  // All relevant faces
	std::vector<double> _cellVolumes;		    // Volumes of the TetGenDelaunay calculated cells

	void ExtractTetGenFaces(const TetGenDelaunay &del);
	void CalculateTetGenVolumes(const TetGenDelaunay &del);
	void OptimizeFace(size_t faceNum); // Optimize the face in place

	//\brief Constructs all the cell objects
	void ConstructCells(const TetGenDelaunay &del);
};

template<typename GhostBusterType>
const double TetGenTessellation<GhostBusterType>::EDGE_RATIO = 1e-5;

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
	ExtractTetGenFaces(del);
	CalculateTetGenVolumes(del);

	for (size_t faceNum = 0; faceNum < _allFaces.size(); faceNum++)
		OptimizeFace(faceNum);

	ConstructCells(del);
}

template<typename GhostBusterType>
void TetGenTessellation<GhostBusterType>::ExtractTetGenFaces(const TetGenDelaunay &del)
{
	const std::vector<Face> _tetgenFaces = del.GetVoronoiFaces();

	_allFaces.assign(_tetgenFaces.size(), boost::none);
	for (size_t cellNum = 0; cellNum < _meshPoints.size(); cellNum++)
	{
		const std::vector<size_t> &faceIndices = del.GetVoronoiCellFaces(cellNum);
		for (std::vector<size_t>::const_iterator it = faceIndices.begin(); it != faceIndices.end(); it++)
			if (!_allFaces[*it].is_initialized())
				_allFaces[*it] = _tetgenFaces[*it];
	}
}

template<typename GhostBusterType>
void TetGenTessellation<GhostBusterType>::CalculateTetGenVolumes(const TetGenDelaunay &del)
{
	_cellVolumes.clear();

	for (size_t cellNum = 0; cellNum < _meshPoints.size(); cellNum++)
	{
		std::vector<const Face *> facePointers;
		const std::vector<size_t> &faceIndices = del.GetVoronoiCellFaces(cellNum);
		for (std::vector<size_t>::const_iterator itFace = faceIndices.begin(); itFace != faceIndices.end(); itFace++)
			facePointers.push_back(&_allFaces[*itFace].value());

		double volume;
		Vector3D CoM;
		CalculateCellDimensions(facePointers, volume, CoM);
		_cellVolumes.push_back(volume);
	}
}

//\brief Optimize a face
//\remark Removes edges that are too short, and then remove faces that are degenerate
//TODO: Remove faces whose area is too small
// We compare each edge to the geometric average of the radius of the two adjacent cells
template<typename GhostBusterType>
void TetGenTessellation<GhostBusterType>::OptimizeFace(size_t faceNum)
{
	if (!_allFaces[faceNum].is_initialized())
		return;

	const Face &face = _allFaces[faceNum].value();
	BOOST_ASSERT(face.NumNeighbors() == 2); // All our faces should have two neighbors, although one of them may be a ghost cell

	int neighbor1 = face.Neighbor1().value().GetCell();
	int neighbor2 = face.Neighbor2().value().GetCell();
	if (neighbor1 > neighbor2)
		std::swap(neighbor1, neighbor2);

	BOOST_ASSERT(neighbor1 < _meshPoints.size()); // At least one neighbor shouldn't be for a ghost
	double volume1 = _cellVolumes[neighbor1];
	double volume2 = neighbor2 < _meshPoints.size() ? _cellVolumes[neighbor2] : volume1;
	double combined = pow(volume1 * volume2, 1.0 / 6.0);  // Cubic Root of geometric average of volumes

	double threshold = combined * EDGE_RATIO;
	double threshold2 = threshold * threshold; // We combine the threshold to the distance, this saves the sqrt operations
	std::vector<VectorRef> vertices; // Vertices we keep in our face

	VectorRef prevVertex = face.vertices.back();
	for (std::vector<VectorRef>::const_iterator it = face.vertices.begin(); it != face.vertices.end(); it++)
	{
		double distance2 = abs2(**it - *prevVertex);
		if (distance2 > threshold2)  // Keep this vertex
		{
			vertices.push_back(*it);
			prevVertex = *it;
		}
	}

	if (vertices.size() == face.vertices.size())  // No change in the face
		return;

	if (vertices.size() < 3)  // Degenerate face
		_allFaces[faceNum] = boost::none;
	else
		_allFaces[faceNum].value().vertices = vertices;
}

template<typename GhostBusterType>
void TetGenTessellation<GhostBusterType>::ConstructCells(const TetGenDelaunay &del)
{
	for (size_t cellNum = 0; cellNum < _meshPoints.size(); cellNum++)
	{
		std::vector<size_t> ourFaceIndices;

		const std::vector<size_t> &tetGenFaceNums = del.GetVoronoiCellFaces(cellNum);
		for (std::vector<size_t>::const_iterator it = tetGenFaceNums.begin(); it != tetGenFaceNums.end(); it++)
		{
			if (!_allFaces[*it].is_initialized())  // Face was removed
				continue;

			const Face &tetGenFace = _allFaces[*it].value();
			size_t ourFaceIndex = _faces.StoreFace(tetGenFace.vertices);
			Face &ourFace = _faces.GetFace(ourFaceIndex);

			ourFaceIndices.push_back(ourFaceIndex);
			ourFace.AddNeighbor(cellNum);
		}

		double volume;
		Vector3D centerOfMass;

		// Convert to face pointers now, because during the previous loop the array may have been reallocated
		// and faces may have been moved in memory
		std::vector<const Face *> ourFacePointers;
		for (std::vector<size_t>::const_iterator it = ourFaceIndices.begin(); it != ourFaceIndices.end(); it++)
			ourFacePointers.push_back(&_faces.GetFace(*it));
		CalculateCellDimensions(ourFacePointers, volume, centerOfMass);

		_cells[cellNum] = Cell(ourFaceIndices, volume, _meshPoints[cellNum], centerOfMass);
		_allCMs[cellNum] = centerOfMass;
	}
}

#endif