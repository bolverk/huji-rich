#include "VoroPlusPlus.hpp"
#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include "../Utilities/assert.hpp"

using namespace std;

VoroPlusPlus::VoroPlusPlus()
{

}

VoroPlusPlus::~VoroPlusPlus()
{

}

void VoroPlusPlus::Initialise(vector<Vector3D> const& points, const OuterBoundary3D &bc)
{
	_boundary = bc;
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
		if (face.Neighbor1()->GetCell() == index)
			neighbors.push_back(face.Neighbor2()->GetCell());
		else if (face.Neighbor2()->GetCell() == index)
			neighbors.push_back(face.Neighbor1()->GetCell());
		else
			BOOST_ASSERT(false); // One of the neighbors must be us!
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


Vector3D VoroPlusPlus::Normal(size_t faceIndex) const
{
	Face face = GetFace(faceIndex);
	Cell cell1 = _cells[face.Neighbor1()->GetCell()];
	Cell cell2 = _cells[face.Neighbor2()->GetCell()];

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

	Face face(vertices);
	index = _faces.size();
	_faces.push_back(face);
	return index;
}

void VoroPlusPlus::FaceStore::Clear()
{
	_faces.clear();
}

void VoroPlusPlus::RunVoronoi()
{
	ClearData();
	auto container = BuildContainer();
	ExtractResults(*container);
}

void VoroPlusPlus::ClearData()
{
	_faces.Clear();
	_cells.clear();
	_cells.resize(_meshPoints.size());
	_allCMs.clear();
	_allCMs.resize(_meshPoints.size());
}

boost::shared_ptr<voro::container> VoroPlusPlus::BuildContainer()
{
	boost::shared_ptr<voro::container> container(new voro::container(
		_boundary.BackLowerLeft().x, _boundary.FrontUpperRight().x,
		_boundary.BackLowerLeft().y, _boundary.FrontUpperRight().y,
		_boundary.BackLowerLeft().z, _boundary.FrontUpperRight().z,
		20, 20, 20, false, false, false, 20));
	int num = 0;
	for (auto it = _meshPoints.begin(); it != _meshPoints.end(); it++)
		container->put(num++, it->x, it->y, it->z);
	return container;
}

void VoroPlusPlus::ExtractResults(voro::container &container)
{
	int cellNum;
	voro::c_loop_all looper(container);

	cellNum = 0;
	looper.start();
	do
	{
		int cellIndex;
		double x, y, z, r;
		looper.pos(cellIndex, x, y, z, r);

		// Create the Cell object
		voro::voronoicell v_cell;
		container.compute_cell(v_cell, looper);
		Cell cell(_meshPoints[cellIndex], v_cell, _faces);
		_cells[cellIndex] = cell;

		// Set the face neighbors
		auto faceIndices = cell.GetFaces();
		for (auto it = faceIndices.begin(); it != faceIndices.end(); it++)
			_faces.GetFace(*it).AddNeighbor(cellIndex); // asserts if fails, no need to add our own assertion

		// And finally, the center of mass
		_allCMs[cellIndex] = cell.GetCenterOfMass();

		cellNum++;
	} while (looper.inc());
	
	BOOST_ASSERT(cellNum == _meshPoints.size());
#ifdef _DEBUG
	for (auto it = _cells.begin(); it != _cells.end(); it++)
		BOOST_ASSERT(!it->empty()); // Make sure we've touched every cell
#endif
}

VoroPlusPlus::Cell::Cell(Vector3D meshPoint, voro::voronoicell &vcell, FaceStore &store)
{
	auto vertices = ExtractAllVertices(meshPoint, vcell);

	vector<int> f_verts;
	vcell.face_vertices(f_verts);
	// f_verts contains a list of faces, each face starts with the number of vertices in the face, followed by
	// the indices of each of the vertices.

	int index = 0;
	for (int faceNum = 0; faceNum < vcell.number_of_faces(); faceNum++)
	{
		int numVertices = f_verts[index++];
		vector<Vector3D> faceVertices;
		for (int i = 0; i < numVertices; i++)
			faceVertices.push_back(vertices[f_verts[index++]]);
		size_t faceIndex = store.StoreFace(faceVertices);
		_faces.push_back(faceIndex);
	}
	BOOST_ASSERT(index == f_verts.size());

	_volume = vcell.volume();
	double cx, cy, cz;
	vcell.centroid(cx, cy, cz);
	_centerOfMass = Vector3D(cx, cy, cz) + meshPoint;
	_center = meshPoint;
	// Volume = (4/3) * Width ^ 3
	_width = pow(3 / 4 * _volume, 1 / 3);
}

vector<Vector3D> VoroPlusPlus::Cell::ExtractAllVertices(Vector3D meshPoint, voro::voronoicell &vcell)
{
	vector<double> voro_vertices;
	vector<Vector3D> our_vertices;

	// Extract the vertices
	vcell.vertices(voro_vertices);
	BOOST_ASSERT(voro_vertices.size() % 3 == 0);

	our_vertices.reserve(voro_vertices.size() / 3);
	for (int i = 0; i < voro_vertices.size(); i += 3)
	{
		Vector3D v(voro_vertices[i], voro_vertices[i + 1], voro_vertices[i + 2]);
		v += meshPoint;
		our_vertices.push_back(v);
	}

	return our_vertices;
}