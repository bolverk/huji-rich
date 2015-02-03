/* \file Delaunay.hpp
   \brief An abstract class that encapsulates the 3D Delaunay calculations
   \Author Itay Zandbank
 */

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
		Vector3D _vec1, _vec2;

	public:
		Edge(Vector3D vec1, Vector3D vec2)
		{
			_vec1 = (vec1 < vec2) ? vec1 : vec2;
			_vec2 = (vec1 < vec2) ? vec2 : vec1;
		}

		const Vector3D &vec1() const { return _vec1; }
		const Vector3D &vec2() const { return _vec2; }
	};

protected:
	std::vector<Vector3D> _points;
	Tetrahedron _bigTetrahedron;

	std::vector<Tetrahedron> _tetrahedra;
	
	typedef std::map<Edge, std::vector<int>> EdgeMap;
	EdgeMap _edges;

	virtual void FillEdges();
	virtual void RunDelaunay() = 0;

public:
	Delaunay(const std::vector<Vector3D> &points, const Tetrahedron &bigTetrahedron);
	void Run();

	bool IsBigTetrahedron(const Vector3D &pt)
	{
		return std::find(_bigTetrahedron.vertices().begin(), 
			_bigTetrahedron.vertices().end(), pt) != _bigTetrahedron.vertices().end();
	}

	const std::vector<Tetrahedron> &Tetrahedra() const { return _tetrahedra; }
	const std::vector<int> &EdgeNeighbors(const Vector3D &vec1, const Vector3D &vec2) const;
};

bool operator==(const Delaunay::Edge &edge1, const Delaunay::Edge &edge2);
bool operator<(const Delaunay::Edge &edge1, const::Delaunay::Edge &edge2);