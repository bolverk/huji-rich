/* \file Delaunay.hpp
   \brief An abstract class that encapsulates the 3D Delaunay calculations
   \Author Itay Zandbank
 */

#ifndef DELAUNAY_HPP
#define DELAUNAY_HPP 1

#include <vector>
#include <list>
#include <map>
#include "../GeometryCommon/Vector3D.hpp"
#include "../GeometryCommon/Tetrahedron.hpp"

class Delaunay
{
public:
	class Edge // This class is public because of the operators required by std::map
	{
	private:
		VectorRef _vec1, _vec2;

	public:
		Edge(VectorRef vec1, VectorRef vec2)
		{
			_vec1 = (vec1 < vec2) ? vec1 : vec2;
			_vec2 = (vec1 < vec2) ? vec2 : vec1;
		}

		const VectorRef &vec1() const { return _vec1; }
		const VectorRef &vec2() const { return _vec2; }
	};

protected:
	std::vector<VectorRef> _points;
	Tetrahedron _bigTetrahedron;

	std::vector<Tetrahedron> _tetrahedra;
	
	typedef std::map<Edge, std::vector<int>> EdgeMap;
	EdgeMap _edges;  // Edge -> list of tetrahedra it appears in

	typedef std::map <VectorRef, std::vector<int>> VertexMap;
	VertexMap _vertices; // Vector->list of tetrahedra it appears in

	virtual void FillVertices();
	virtual void FillEdges();
	virtual void RunDelaunay() = 0;

public:
	Delaunay(const std::vector<VectorRef> &points, const Tetrahedron &bigTetrahedron);
	void Run();

	bool IsBigTetrahedron(const VectorRef &pt) const
	{
		return std::find(_bigTetrahedron.vertices().begin(), 
			_bigTetrahedron.vertices().end(), pt) != _bigTetrahedron.vertices().end();
	}

	//\brief A vector of all the tetrahedra
	const std::vector<Tetrahedron> &Tetrahedra() const { return _tetrahedra; }

	//\brief Returns the list of tetrahedra that touch the edge 
	const std::vector<int> &EdgeNeighbors(const VectorRef vec1, const VectorRef vec2) const;
	
	//\brief Returns the list of tetrahedra that touch the vertex
	const std::vector<int> &VertexNeighbors(const VectorRef v) const;

	size_t NumTetrahedra() const { return _tetrahedra.size(); }
	const Tetrahedron& operator[](size_t index) const { return _tetrahedra[index]; }

	const Tetrahedron& BigTetrahedron() const { return _bigTetrahedron; }
	const std::vector<VectorRef> &InputPoints() const { return _points; }
};

bool operator==(const Delaunay::Edge &edge1, const Delaunay::Edge &edge2);
bool operator<(const Delaunay::Edge &edge1, const::Delaunay::Edge &edge2);

#endif // DELAUNAY_HPP