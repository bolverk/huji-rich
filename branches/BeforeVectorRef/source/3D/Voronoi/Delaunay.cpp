/*  \file Delaunay.cpp
\brief An abstract class that encapsulates the 3D Delaunay calculations
\Author Itay Zandbank
*/

#include "../Utilities/assert.hpp"
#include "Delaunay.hpp"

Delaunay::Delaunay(const std::vector<Vector3D> &points, const Tetrahedron &bigTetrahedron) 
	: _points(points), _bigTetrahedron(bigTetrahedron)
{

}

/* Fills the _edges map after the Tetrahedra have been set */
void Delaunay::FillEdges()
{
	_edges.clear();

	for (int i = 0; i < _tetrahedra.size(); i++)
	{
		const Tetrahedron &t = _tetrahedra[i];
		for (int iv = 0; iv < 3; iv++)  // Go over all 6 edges
			for (int jv = iv + 1; jv < 4; jv++)
			{
				Edge edge(t[iv], t[jv]);
				EdgeMap::iterator existing = _edges.find(edge);
				if (existing == _edges.end())
					_edges[edge] = vector<int>();

				_edges[edge].push_back(i);
			}
	}
}

void Delaunay::FillVertices()
{
	_vertices.clear();

	for (int i = 0; i < _tetrahedra.size(); i++)
		for (int j = 0; j < 4; j++)
		{
			Vector3D v = _tetrahedra[i][j];
			VertexMap::iterator existing = _vertices.find(v);
			if (existing == _vertices.end())
				_vertices[v] = vector<int>();

			_vertices[v].push_back(i);
		}
}

const std::vector<int> &Delaunay::EdgeNeighbors(const Vector3D &vec1, const Vector3D &vec2) const
{
	Edge edge(vec1, vec2);
	EdgeMap::const_iterator it = _edges.find(edge);
	BOOST_ASSERT(it != _edges.end());

	return it->second;
}

const std::vector<int> &Delaunay::VertexNeighbors(const Vector3D &v) const
{
	VertexMap::const_iterator it = _vertices.find(v);
	BOOST_ASSERT(it != _vertices.end());

	return it->second;
}

bool operator==(const Delaunay::Edge &edge1, const Delaunay::Edge &edge2)
{
	return edge1.vec1() == edge2.vec1() && edge1.vec2() == edge2.vec2();
}

bool operator<(const Delaunay::Edge &edge1, const::Delaunay::Edge &edge2)
{
	if (edge1.vec1() < edge2.vec1())
		return true;
	else if (edge1.vec1() == edge2.vec1())
		return edge1.vec2() < edge2.vec2();
	else
		return false;
}

void Delaunay::Run()
{
	_tetrahedra.clear();
	_edges.clear();
	_vertices.clear();
	RunDelaunay();
	BOOST_ASSERT(_edges.size() == 0 && _vertices.size() == 0); // You're not supposed to fill those
	FillVertices();
	FillEdges();
}