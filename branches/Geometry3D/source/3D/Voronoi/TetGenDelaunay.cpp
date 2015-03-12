/*
\file TetGenDelaunay.cpp
\brief a TetGen based implementation of the Delaunay abstract class
\author Itay Zandbank
*/

#include "TetGenDelaunay.hpp"
#include <tetgen.h>
#include <vector>
#include <unordered_map>
#include <algorithm>

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

	const VectorRef GetPoint(size_t offset);

public:
	TetGenImpl(TetGenDelaunay &delaunay) : _delaunay(delaunay) { }
	void Run();
};

TetGenDelaunay::TetGenDelaunay(const std::vector<VectorRef> &points, const Tetrahedron &bigTetrahedron, bool runVoronoi)
	: Delaunay(points, bigTetrahedron), _runVoronoi(runVoronoi)
{
}

void TetGenDelaunay::RunDelaunay()
{
	TetGenImpl impl(*this);
	impl.Run();
}

const VectorRef TetGenImpl::GetPoint(size_t offset)
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
	for (vector<VectorRef>::const_iterator it = _delaunay._points.begin(); it != _delaunay._points.end(); it++)
	{
		VectorRef vec = *it;
		in.pointlist[offset++] = vec->x;
		in.pointlist[offset++] = vec->y;
		in.pointlist[offset++] = vec->z;
	}

	// Copy the big tetrahedron
	for (int i = 0; i < 4; i++)
	{
		in.pointlist[offset++] = _delaunay._bigTetrahedron[i]->x;
		in.pointlist[offset++] = _delaunay._bigTetrahedron[i]->y;
		in.pointlist[offset++] = _delaunay._bigTetrahedron[i]->z;
	}
}

void TetGenImpl::CallTetGen()
{
	tetgenbehavior b;
	std::string commandLine = "nQ";
	if (_delaunay._runVoronoi)
		commandLine += "v";
	b.parse_commandline((char *)commandLine.c_str()); // parse_command line expects a non-const pointer. Blach.
	tetrahedralize(&b, &in, &out);
}

void TetGenImpl::CopyResults()
{
	_delaunay._tetrahedraNeighbors.clear();

	int offset = 0;
	for (int i = 0; i < out.numberoftetrahedra; i++)
	{
		vector<VectorRef> vertices;
		vector<size_t> neighbors;
		for (int j = 0; j < 4; j++)
		{
			int ptIndex = out.tetrahedronlist[offset];
			vertices.push_back(GetPoint(ptIndex));

			int neighbor = out.neighborlist[offset];
			if (neighbor >= 0)
				neighbors.push_back(out.neighborlist[offset]);
			offset++;
		}
		Tetrahedron t(vertices);
		_delaunay._tetrahedra.push_back(t);
		_delaunay._tetrahedraNeighbors.push_back(neighbors);
	}

	if (_delaunay._runVoronoi) // Copy the Voronoi output
	{
		_delaunay._voronoiCellFaces.clear();
		_delaunay._voronoiFaceEdges.clear();
		_delaunay._voronoiVertices.clear();

		for (size_t cellNum = 0; cellNum < out.numberofvcells; cellNum++)
		{
			int faceCount = out.vcelllist[cellNum][0];
			vector<int> faces;
			faces.assign(out.vcelllist[cellNum] + 1, out.vcelllist[cellNum] + faceCount + 1);
			_delaunay._voronoiCellFaces.push_back(faces);
		}

		for (size_t faceNum = 0; faceNum < out.numberofvfacets; faceNum++)
		{
			const tetgenio::vorofacet &face = out.vfacetlist[faceNum];
			int edgeCount = face.elist[0];
			vector<int> edges;
			edges.assign(face.elist + 1, face.elist + edgeCount + 1);
			_delaunay._voronoiFaceEdges.push_back(edges);

			_delaunay._voronoiFaceNeighbors.push_back(pair<int, int>(face.c1, face.c2));
		}

		for (size_t edgeNum = 0; edgeNum < out.numberofvedges; edgeNum++)
		{
			pair<int, int> edge(out.vedgelist[edgeNum].v1, out.vedgelist[edgeNum].v2);
			_delaunay._voronoiEdges.push_back(edge);
		}

		double *ptr = out.vpointlist;
		for (size_t vertexNum = 0; vertexNum < out.numberofvpoints; vertexNum++)
		{
			Vector3D vertex(ptr[0], ptr[1], ptr[2]);
			_delaunay._voronoiVertices.push_back(vertex);
			ptr += 3;
		}
	}
}

/* Fills the _edges map after the Tetrahedra have been set */
void TetGenDelaunay::FillEdges()
{
	_edges.clear();

	// First step, fill the edge map with a list of tetrahedra that touch the edge, with no particular order
	for (size_t i = 0; i < _tetrahedra.size(); i++)
	{
		const Tetrahedron &t = _tetrahedra[i];
		for (int iv = 0; iv < 3; iv++)  // Go over all 6 edges
			for (int jv = iv + 1; jv < 4; jv++)
			{
				Edge edge(t[iv], t[jv]);
				EdgeMap::iterator existing = _edges.find(edge);
				if (existing == _edges.end())
					_edges[edge] = vector<size_t>();

				_edges[edge].push_back(i);
			}
	}

	// Second step, order the edge map properly
	for (EdgeMap::iterator it = _edges.begin(); it != _edges.end(); it++)
	{
		vector<size_t> tetrahedra = it->second;
		vector<size_t> ordered = OrderNeighbors(tetrahedra);
		it->second = ordered;
	}
}

/* Fills the neighbors. This actually does nothing, because CopyResults already does this */
void TetGenDelaunay::FillNeighbors()
{
	// Already done by TetGenImpl::CopyResults
}

// Order the neighbors of the edge so that they touch each-other
// If there edge is not a boundary edge, it's easy since there tetrahedra form a cycle
//
// Boundary edges, however, don't have a cycle of neighbors. They tetrahedra do form a path,
// So if we find out we can't go any further to one side, we just start walking to the other side.
// This is done by reversing the neighbors found so far and continuing.

vector<size_t> TetGenDelaunay::OrderNeighbors(const vector<size_t> &tetrahedra)
{
	BOOST_ASSERT(tetrahedra.size() > 0);
	vector<size_t> ordered;
	ordered.reserve(tetrahedra.size());

	bool changedDirection = false; // True if we change directions - so we don't change directions twice
	ordered.push_back(tetrahedra[0]); // Arbitrary starting point

	int index = 1;
	while (index < tetrahedra.size())
	{ 
		// Find the next neighbor
		size_t previous = ordered[index - 1];
		const vector<size_t> &previousNeighbors = TetrahedraNeighbors(previous);
		vector<size_t>::const_iterator it;
		for (it = previousNeighbors.begin(); it != previousNeighbors.end(); it++)
		{
			// TODO: This is all O(N^2) - improve if this is too lengthy in the real world
			// See if this neighbor is part of the edge
			if (find(tetrahedra.begin(), tetrahedra.end(), *it) == tetrahedra.end())
				continue; // Nope

			// See if the neighbor was already used
			if (find(ordered.begin(), ordered.end(), *it) == ordered.end())
				break;
		}

		// Found a neighbor, use it
		if (it != previousNeighbors.end())
		{
			ordered.push_back(*it);
			index++; // Move to the next one
		}
		else if (!changedDirection) // Change direction
		{
			reverse(ordered.begin(), ordered.end());
			changedDirection = true;
		}
		else
		{
			BOOST_ASSERT(false); // No neighbor even though we changed direction, this is a bug!
			return vector<size_t>();
		}
	}

	return ordered;
}


vector<size_t> TetGenDelaunay::GetVoronoiCellFaces(size_t cellNum) const
{
	vector<size_t> faces(_voronoiCellFaces[cellNum].begin(), _voronoiCellFaces[cellNum].end()); // Convert to size_t

	return faces;
}

// Construct an ordered list of vertices of the face.
// The face consists of edges, each edge consists of two vertices. Each edge shares
// one vertex with the edge before it. So each edge, we take the vertex that does not appear
// in the edge before it and use it.
//
// An edge list for example:   (1,5)   (5, 8)  (4, 8)  (3 4)   (3 10) (1 10)
// Means the edge consists of vertices 1 5 8 4 3 10
// So each edge translates to one vertex - the one that does not appear in the previous edge
vector<Vector3D> TetGenDelaunay::GetVoronoiFace(size_t faceNum) const
{
	const vector<int> &edgeList = _voronoiFaceEdges[faceNum];

	BOOST_ASSERT(edgeList.back() >= 0); // We shouldn't be asking for an open face!
	if (edgeList.back() == -1)
		return vector<Vector3D>();   // Return an empty face

	vector<int> indices;
	pair<int, int> previous = _voronoiEdges[edgeList.back()];
	for (vector<int>::const_iterator it = edgeList.begin(); it != edgeList.end(); it++)
	{
		pair<int, int> current = _voronoiEdges[*it];

		if (current.first == previous.first || current.first == previous.second)  // current.first appears in the previous edge
			indices.push_back(current.second);
		else if (current.second == previous.first || current.second == previous.second) // current.second appears in the previous edge
			indices.push_back(current.first);
		else
			BOOST_ASSERT(false); // Non-consecutive edges!

		previous = current;
	}

	// Now, convert the list of indices to a list of vectors
	vector<Vector3D> face;
	for (vector<int>::const_iterator it = indices.begin(); it != indices.end(); it++)
	{
		BOOST_ASSERT(*it >= 0); // We shouldn't be asking for any rays!
		if (*it == -1)
			return vector<Vector3D>();
		face.push_back(_voronoiVertices[*it]);
	}

	return face;
}