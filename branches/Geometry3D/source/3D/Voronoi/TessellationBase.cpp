#include "TessellationBase.hpp"
#include <set>

/*
* The FaceStore
*/

bool TessellationBase::FaceStore::FindFace(const vector<VectorRef> &vertices, size_t &index) const
{
	for (index = 0; index < _faces.size(); index++)
		if (_faces[index].IdenticalTo(vertices))
			return true;

	return false;
}

size_t TessellationBase::FaceStore::StoreFace(const vector<VectorRef> &vertices)
{
	size_t index;
	bool exists = FindFace(vertices, index);
	if (exists)
		return index;

	Face face(vertices);
	index = _faces.size();
	_faces.push_back(face);
	return index;
}

void TessellationBase::FaceStore::Clear()
{
	_faces.clear();
}

TessellationBase::Cell::Cell(std::vector<size_t> faces, double volume, double width, Vector3D center, Vector3D centerOfMass) :
	_faces(faces), _volume(volume), _width(width), _center(center), _centerOfMass(centerOfMass)
{

}

void TessellationBase::ClearData()
{
	_faces.Clear();
	_cells.clear();
	_cells.resize(_meshPoints.size());
	_allCMs.clear();
	_allCMs.resize(_meshPoints.size());
	FillPointIndices();
}

void TessellationBase::Initialise(vector<Vector3D> const& points, const OuterBoundary3D &bc)
{
	_boundary = &bc;
	Update(points);
}

size_t TessellationBase::GetPointNo(void) const
{
	return _meshPoints.size();
}

Vector3D TessellationBase::GetMeshPoint(size_t index) const
{
	return _meshPoints[index];
}

Vector3D const& TessellationBase::GetCellCM(size_t index) const
{
	return _cells[index].GetCenterOfMass();
}

/*! \brief Returns the total number of faces
\return Total number of faces
*/
size_t TessellationBase::GetTotalFacesNumber(void) const
{
	return _faces.NumFaces();
}

const Face& TessellationBase::GetFace(size_t index) const
{
	return _faces.GetFace(index);
}

double TessellationBase::GetWidth(size_t index) const
{
	return _cells[index].GetWidth();
}

double TessellationBase::GetVolume(size_t index) const
{
	return _cells[index].GetVolume();
}

vector<size_t>const& TessellationBase::GetCellFaces(size_t index) const
{
	return _cells[index].GetFaces();
}

vector<Vector3D>& TessellationBase::GetMeshPoints(void)
{
	return _meshPoints;
}

vector<size_t> TessellationBase::GetNeighbors(size_t index) const
{
	const Cell& cell = _cells[index];
	vector<size_t> neighbors;

	auto faceIndices = _cells[index].GetFaces();
	for (size_t i = 0; i < faceIndices.size(); i++)
	{
		auto face = GetFace(i);
		if (face.Neighbor1()->GetCell() == index)
			neighbors.push_back(face.Neighbor2()->GetCell());
		else if (face.Neighbor2()->GetCell() == index)
			neighbors.push_back(face.Neighbor1()->GetCell());
		else
			BOOST_ASSERT(false); // One of the neighbors must be us!
	}

	return neighbors;
}


vector<Vector3D>& TessellationBase::GetAllCM()
{
	return _allCMs;
}

void TessellationBase::GetNeighborNeighbors(vector<size_t> &result, size_t point) const
{
	set<size_t> allNeighbors;
	auto neighbors = GetNeighbors(point);
	for (auto it = neighbors.begin(); it != neighbors.end(); it++)
	{
		allNeighbors.insert(*it);
		auto neighbors2 = GetNeighbors(*it);
		allNeighbors.insert(neighbors2.begin(), neighbors2.end());
	}

	// See here: http://stackoverflow.com/a/5034274/871910
	result.clear();
	copy(allNeighbors.begin(), allNeighbors.end(), back_inserter(result));
}


Vector3D TessellationBase::Normal(size_t faceIndex) const
{
	Face face = GetFace(faceIndex);
	Cell cell1 = _cells[face.Neighbor1()->GetCell()];
	Cell cell2 = _cells[face.Neighbor2()->GetCell()];

	return cell1.GetCenterOfMass() - cell2.GetCenterOfMass();
}

/*!
\brief Calculates the velocity of a face
\param p0 The index of the first neighbor
\param p1 The index of the second neighbor
\param v0 The velocity of the first neighbor
\param v1 The velocity of the second neighbor
\return The velocity of the face
*/
Vector3D TessellationBase::CalcFaceVelocity(size_t p0, size_t p1, Vector3D const& v0,
	Vector3D const& v1) const
{
	return Vector3D(); // TODO: Calculate the velocity properly
}

//\brief Check if a cell touches the boundary
//\remarks The cell touches a boundary iff it has a boundary face
bool TessellationBase::NearBoundary(size_t index) const
{
	const vector<size_t> faces = GetCellFaces(index);
	for (vector<size_t>::const_iterator it = faces.begin(); it != faces.end(); it++)
		if (BoundaryFace(*it))
			return true;
	return false;
}

//\brief Checks if a face is a Boundary face
//\remarks a Face is a Boundary face iff it only has one neighbor
bool TessellationBase::BoundaryFace(size_t index) const
{
	Face face = GetFace(index);
	return face.NumNeighbors() == 1;
}

//\brief Fills the pointIndices map
void TessellationBase::FillPointIndices()
{
	_pointIndices.clear();
	for (size_t i = 0; i < _meshPoints.size(); i++)
		_pointIndices[_meshPoints[i]] = i;
}

boost::optional<size_t> TessellationBase::GetPointIndex(const VectorRef pt) const
{
	unordered_map<VectorRef, size_t>::const_iterator it = _pointIndices.find(pt);
	if (it != _pointIndices.end())
		return it->second;

	return boost::none;
}

void TessellationBase::GetTetrahedronIndices(const Tetrahedron &t, boost::optional<size_t> *cells) const
{
	for (int i = 0; i < 4; i++)
		cells[i] = GetPointIndex(t[i]);
}

