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
	return _faces.size();
}

Face const& VoroPlusPlus::GetFace(size_t index) const
{
	return _faces[index];
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
	for (int i = 0; i < faceIndices.size(); i++)
	{
		auto face = _faces[i];
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
	vector<vector<size_t>> v;
	return v;
}

vector<vector<size_t> >const& VoroPlusPlus::GetDuplicatedPoints(void) const
{
	vector<vector<size_t>> v;
	return v;
}

size_t VoroPlusPlus::GetTotalPointNumber(void)const
{
	return _meshPoints.size();
}

vector<Vector3D>& VoroPlusPlus::GetAllCM()
{
	vector<Vector3D> vec(_cells.size());
	for (int i = 0; i < vec.size(); i++)
		vec[i] = _cells[i].GetCenterOfMass();

	return vec;
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


Vector3D VoroPlusPlus::Normal(size_t faceindex) const
{
	Face face = _faces[faceindex];
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
