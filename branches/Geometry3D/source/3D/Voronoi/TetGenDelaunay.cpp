/*
\file TetGenDelaunay.cpp
\brief a TetGen based implementation of the Delaunay abstract class
\author Itay Zandbank
*/

#include "TetGenDelaunay.hpp"
#include <tetgen.h>

using namespace std; 

// This class exists to hide the TetGen library from the callers - all the TetGen dependant code is in this source
// file. Library users will not need to add tetgen.h to the include path, they only need to link against the TetGen library.
class TetGenImpl
{
private:
	TetGenDelaunay &_delaunay;

	tetgenio in, out;

	void InitInput();
	void CallTetGen();
	void CopyResults();

	const Vector3D &GetPoint(size_t offset);

public:
	TetGenImpl(TetGenDelaunay &delaunay) : _delaunay(delaunay) { }
	void Run();
};

TetGenDelaunay::TetGenDelaunay(const std::vector<Vector3D> &points, const Tetrahedron &bigTetrahedron)
	: Delaunay(points, bigTetrahedron)
{
}

void TetGenDelaunay::RunDelaunay()
{
	TetGenImpl impl(*this);
	impl.Run();
}

const Vector3D &TetGenImpl::GetPoint(size_t offset)
{
	BOOST_ASSERT(offset >= 0);
	BOOST_ASSERT(offset < _delaunay._points.size() + 4);

	if (offset < _delaunay._points.size())
		return _delaunay._points[offset];
	offset -= _delaunay._points.size();
	return _delaunay._bigTetrahedron[offset];
}

void TetGenImpl::Run()
{
	InitInput();
	CallTetGen();
	CopyResults();
}

void TetGenImpl::InitInput()
{
	in.firstnumber = 0;
	in.numberofpoints = (int)(4 + _delaunay._points.size());  // 4 for the Big Tetrahedron
	in.pointlist = new REAL[3 * in.numberofpoints];

	int offset = 0;
	// Copy the points
	for (vector<Vector3D>::const_iterator it = _delaunay._points.begin(); it != _delaunay._points.end(); it++)
	{
		in.pointlist[offset++] = it->x;
		in.pointlist[offset++] = it->y;
		in.pointlist[offset++] = it->z;
	}

	// Copy the big tetrahedron
	for (int i = 0; i < 4; i++)
	{
		in.pointlist[offset++] = _delaunay._bigTetrahedron[i].x;
		in.pointlist[offset++] = _delaunay._bigTetrahedron[i].y;
		in.pointlist[offset++] = _delaunay._bigTetrahedron[i].z;
	}
}

void TetGenImpl::CallTetGen()
{
	tetgenbehavior b;
	b.parse_commandline("efnvQ");
	tetrahedralize(&b, &in, &out);
}

void TetGenImpl::CopyResults()
{
	int offset = 0;
	for (int i = 0; i < out.numberoftetrahedra; i++)
	{
		vector<Vector3D> vertices;
		vector<int> neighbors;
		for (int j = 0; j < 4; j++)
		{
			int ptIndex = out.tetrahedronlist[offset];
			vertices.push_back(GetPoint(ptIndex));

			neighbors.push_back(out.neighborlist[offset]);
			offset++;
		}
		Tetrahedron t(vertices);
		_delaunay._tetrahedra.push_back(t);
	}
}
