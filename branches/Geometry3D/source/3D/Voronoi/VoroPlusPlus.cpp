#include "VoroPlusPlus.hpp"
#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include "../Utilities/assert.hpp"
#include <voro++.hh>

using namespace std;

// This class contains all the voro++ references. This means callers to this library don't need to include voro++.hh,
// they just need to link against the Voro++ library.
class VoroPlusPlusImpl
{
private:
	VoroPlusPlus &_voro;

public:
	VoroPlusPlusImpl(VoroPlusPlus &voro);
	void Run();

private:
	int FindMeshPoint(voro::c_loop_all looper);
	boost::shared_ptr<voro::container> BuildContainer();
	void ExtractResults(voro::container &container);

	TessellationBase::Cell CreateCell(Vector3D meshPoint, voro::voronoicell &vcell);
	vector<Vector3D> ExtractAllVertices(Vector3D meshPoint, voro::voronoicell &vcell);
};

VoroPlusPlus::VoroPlusPlus()
{

}

VoroPlusPlus::~VoroPlusPlus()
{

}

void VoroPlusPlus::Update(vector<Vector3D> const& points)
{
	_meshPoints = points;
	RunVoronoi();
}


Tessellation3D* VoroPlusPlus::clone(void) const
{
	return new VoroPlusPlus(*this);
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


bool VoroPlusPlus::IsGhostPoint(size_t index) const
{
	// TODO: Treat ghosts
	return false;
}

void VoroPlusPlus::RunVoronoi()
{
	ClearData();
	VoroPlusPlusImpl impl(*this);
	impl.Run();
}

VoroPlusPlusImpl::VoroPlusPlusImpl(VoroPlusPlus &voro) : _voro(voro) { }

void VoroPlusPlusImpl::Run()
{
	auto container = BuildContainer();
	ExtractResults(*container);
}

boost::shared_ptr<voro::container> VoroPlusPlusImpl::BuildContainer()
{
	boost::shared_ptr<voro::container> container(new voro::container(
		_voro._boundary->BackLowerLeft().x, _voro._boundary->FrontUpperRight().x,
		_voro._boundary->BackLowerLeft().y, _voro._boundary->FrontUpperRight().y,
		_voro._boundary->BackLowerLeft().z, _voro._boundary->FrontUpperRight().z,
		20, 20, 20, false, false, false, 20));
	int num = 0;
	for (auto it = _voro._meshPoints.begin(); it != _voro._meshPoints.end(); it++)
		container->put(num++, it->x, it->y, it->z);
	return container;
}

void VoroPlusPlusImpl::ExtractResults(voro::container &container)
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
		// Cell cell(_meshPoints[cellIndex], v_cell, _faces);
		VoroPlusPlus::Cell cell = CreateCell(_voro._meshPoints[cellIndex], v_cell);
		_voro._cells[cellIndex] = cell;

		// Set the face neighbors
		auto faceIndices = cell.GetFaces();
		for (auto it = faceIndices.begin(); it != faceIndices.end(); it++)
			_voro._faces.GetFace(*it).AddNeighbor(cellIndex); // asserts if fails, no need to add our own assertion

		// And finally, the center of mass
		_voro._allCMs[cellIndex] = cell.GetCenterOfMass();

		cellNum++;
	} while (looper.inc());
	
	BOOST_ASSERT(cellNum == _voro._meshPoints.size());
#ifdef _DEBUG
	for (auto it = _voro._cells.begin(); it != _voro._cells.end(); it++)
		BOOST_ASSERT(!it->empty()); // Make sure we've touched every cell
#endif
}

TessellationBase::Cell VoroPlusPlusImpl::CreateCell(Vector3D meshPoint, voro::voronoicell &vcell)
{
	auto vertices = ExtractAllVertices(meshPoint, vcell);

	vector<int> f_verts;
	vcell.face_vertices(f_verts);
	// f_verts contains a list of faces, each face starts with the number of vertices in the face, followed by
	// the indices of each of the vertices.

	int index = 0;
	vector<size_t> faces;
	for (int faceNum = 0; faceNum < vcell.number_of_faces(); faceNum++)
	{
		int numVertices = f_verts[index++];
		vector<VectorRef> faceVertices;
		for (int i = 0; i < numVertices; i++)
			faceVertices.push_back(vertices[f_verts[index++]]);
		size_t faceIndex = _voro._faces.StoreFace(faceVertices);
		faces.push_back(faceIndex);
	}
	BOOST_ASSERT(index == f_verts.size());

	double volume = vcell.volume();
	double cx, cy, cz;
	vcell.centroid(cx, cy, cz);
	auto centerOfMass = Vector3D(cx, cy, cz) + meshPoint;
	//_center = meshPoint;
	// Volume = (4/3) * Width ^ 3
	double width = pow(3 / 4 * volume, 1 / 3);

	return TessellationBase::Cell(faces, volume, width, meshPoint, centerOfMass);
}

vector<Vector3D> VoroPlusPlusImpl::ExtractAllVertices(Vector3D meshPoint, voro::voronoicell &vcell)
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
