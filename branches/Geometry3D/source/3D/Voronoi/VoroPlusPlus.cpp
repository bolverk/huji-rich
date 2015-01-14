#include "VoroPlusPlus.hpp"
#include <vector>
#include <set>

using namespace std;

VoroPlusPlus::VoroPlusPlus()
{

}

VoroPlusPlus::~VoroPlusPlus()
{

}

void VoroPlusPlus::Initialise(vector<Vector3D> const& points, OuterBoundary3D const* bc)
{
	// TODO: Store boundry conditions
	Update(points);
}

void VoroPlusPlus::Update(vector<Vector3D> const& points)
{
	_meshPoints = points;
	RunVoronoi();
}

size_t VoroPlusPlus::GetPointNo(void) const
{
	return _meshPoints.size();
}

Vector3D VoroPlusPlus::GetMeshPoint(size_t index) const
{
	return _meshPoints[index];
}

Vector3D const& VoroPlusPlus::GetCellCM(size_t index) const
{
	return _cells[index].GetCenterOfMass();
}

/*! \brief Returns the total number of faces
\return Total number of faces
*/
size_t VoroPlusPlus::GetTotalFacesNumber(void) const
{
	return _faces.NumFaces();
}

const Face& VoroPlusPlus::GetFace(size_t index) const
{
	return _faces.GetFace(index);
}

double VoroPlusPlus::GetWidth(size_t index) const
{
	return _cells[index].GetWidth();
}

double VoroPlusPlus::GetVolume(size_t index) const
{
	return _cells[index].GetVolume();
}

vector<size_t>const& VoroPlusPlus::GetCellFaces(size_t index) const
{
	return _cells[index].GetFaces();
}

vector<Vector3D>& VoroPlusPlus::GetMeshPoints(void)
{
	return _meshPoints;
}

vector<size_t> VoroPlusPlus::GetNeighbors(size_t index) const
{
	const Cell& cell = _cells[index];
	vector<size_t> neighbors;

	auto faceIndices = _cells[index].GetFaces();
	for (size_t i = 0; i < faceIndices.size(); i++)
	{
		auto face = GetFace(i);
		if (face.neighbors.first == index)
			neighbors.push_back(face.neighbors.second);
		else if (face.neighbors.second == index)
			neighbors.push_back(face.neighbors.first);
		else
			assert(false); // One of the neighbors must be us!
	}

	return neighbors;
}

Tessellation3D* VoroPlusPlus::clone(void) const
{
	return new VoroPlusPlus(*this);
}

bool VoroPlusPlus::NearBoundary(size_t index) const
{
	// TODO: This
	return false;
}

bool VoroPlusPlus::BoundaryFace(size_t index) const
{
	// TODO: This
	return false;
}

vector<vector<size_t> >& VoroPlusPlus::GetDuplicatedPoints()
{
	static vector<vector<size_t>> v;
	return v;
}

vector<vector<size_t> >const& VoroPlusPlus::GetDuplicatedPoints(void) const
{
	static vector<vector<size_t>> v;
	return v;
}

size_t VoroPlusPlus::GetTotalPointNumber(void)const
{
	return _meshPoints.size();
}

vector<Vector3D>& VoroPlusPlus::GetAllCM()
{
	return _allCMs;
}

void VoroPlusPlus::GetNeighborNeighbors(vector<size_t> &result, size_t point) const
{
	set<int> allNeighbors;
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


Vector3D VoroPlusPlus::Normal(size_t faceIndex) const
{
	Face face = GetFace(faceIndex);
	Cell cell1 = _cells[face.neighbors.first];
	Cell cell2 = _cells[face.neighbors.second];

	return cell1.GetCenterOfMass() - cell2.GetCenterOfMass();
}

bool VoroPlusPlus::IsGhostPoint(size_t index) const
{
	// TODO: Treat ghosts
	return false;
}

/*!
\brief Calculates the velocity of a face
\param p0 The index of the first neighbor
\param p1 The index of the second neighbor
\param v0 The velocity of the first neighbor
\param v1 The velocity of the second neighbor
\return The velocity of the face
*/
Vector3D VoroPlusPlus::CalcFaceVelocity(size_t p0, size_t p1, Vector3D const& v0,
	Vector3D const& v1) const
{
	return Vector3D(); // TODO: Calculate the velocity properly
}

/*
 * The FaceStore
 */

bool VoroPlusPlus::FaceStore::FindFace(const vector<Vector3D> &vertices, size_t &index) const
{
	for (index = 0; index < _faces.size(); index++)
		if (_faces[index].IdenticalTo(vertices))
			return true;

	return false;
}

size_t VoroPlusPlus::FaceStore::StoreFace(const vector<Vector3D> &vertices)
{
	size_t index;
	bool exists = FindFace(vertices, index);
	if (exists)
		return index;

	Face face(vertices, (size_t)-1, (size_t)-1);
	index = _faces.size();
	_faces.push_back(face);
	return index;
}

void VoroPlusPlus::RunVoronoi()
{
	auto container = BuildContainer();
	ExtractResults(container);
}

voro::container VoroPlusPlus::BuildContainer()
{
	Vector3D min, max;

	min = _meshPoints[0];
	max = _meshPoints[0];

	// Find the boundry
	for (auto it = _meshPoints.begin(); it != _meshPoints.end(); it++)
	{
		min.x = std::min(min.x, it->x);
		min.y = std::min(min.y, it->y);
		min.z = std::min(min.z, it->z);
		max.x = std::max(max.x, it->x);
		max.y = std::max(max.y, it->y);
		max.z = std::max(max.z, it->z);
	}

	voro::container container(min.x, max.x, min.y, max.y, min.z, max.z, 20, 20, 20, false, false, false, 20);
	int num = 0;
	for (auto it = _meshPoints.begin(); it != _meshPoints.end(); it++)
		container.put(num++, it->x, it->y, it->z);

	return container;
}

void VoroPlusPlus::ExtractResults(voro::container container)
{
	// TODO: Extract Results
}