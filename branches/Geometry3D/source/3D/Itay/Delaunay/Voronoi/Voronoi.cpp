// Voronoi.cpp : Defines the entry point for the console application.
//

#include "tetgen.h"
#include "GeometryCommon/Mat44.hpp"
#include "GeometryCommon/Vector3D.hpp"
#include "Voronoi/VoroPlusPlus.hpp"
#include "Voronoi/TetGenDelaunay.hpp"
#include <vector>
#include <iostream>
#include <set>
#include <unordered_set>
#include <algorithm>
#include "GeometryCommon/Tetrahedron.hpp"
#include <map>
#include <sstream>
#include <string>
#include "GeometryCommon/Subcube.hpp"
#include "Voronoi/DelaunayVoronoi.hpp"

using namespace std;

class VectorNamer
{
private:
	map<Vector3D, string> _vectors;
	int _last;

public:
	VectorNamer()
	{
		_last = 0;
	}

	string GetName(Vector3D vec, string prefix="V")
	{
		auto inMap = _vectors.find(vec);
		if (inMap != _vectors.end())
			return inMap->second;

		stringstream strm;
		strm << prefix << _last++;
		_vectors[vec] = strm.str();
		return strm.str();
	}

	map<Vector3D, string>::const_iterator begin() const { return _vectors.cbegin(); }
	map<Vector3D, string>::const_iterator end() const { return _vectors.cend(); }
};

set<Vector3D> UseVoroPlusPlus(const vector<Vector3D>& points, const OuterBoundary3D &boundary);
set<Vector3D> UseTetGen(const vector<Vector3D>& points, const OuterBoundary3D &boundary);

Tetrahedron FindBigTetrahedron(const OuterBoundary3D &boundary);
unordered_set<int> FindOuterTetrahedra(const Delaunay &del);
unordered_set<int> FindEdgeTetrahedra(const Delaunay &del, const unordered_set<int>& outerTetrahedra);
map<Vector3D, unordered_set<Subcube>> FindHullBreaches(const Delaunay &del, const unordered_set<int>& edgeTetrahedra, const unordered_set<int> &outerTetrahedra, const OuterBoundary3D &boundary);
set<Vector3D> GetGhostPoints(const Delaunay &del, const OuterBoundary3D &boundary);
set<Vector3D> BruteForceGhosts(const Delaunay &del, const OuterBoundary3D &boundary);

template<typename T>
bool contains(const unordered_set<T> &set, const T &val)
{
	return set.find(val) != set.end();
}

VectorNamer namer;

int main()
{
	vector<Vector3D> vertices{ Vector3D(90, 34, 89),
		Vector3D(21, 3, 78),
		Vector3D(76, 35, 74),
		Vector3D(28, 4, 7),
		Vector3D(65, 60, 22),
		Vector3D(59, 92, 5),
		Vector3D(84, 0, 32) };
	RectangularBoundary3D boundary(Vector3D(200, 200, 200), Vector3D(-200, -200, -200));

	for (auto vec : vertices)
		namer.GetName(vec, "C");

	auto voroVertices = UseVoroPlusPlus(vertices, boundary);
	auto delaunayVertices = UseTetGen(vertices, boundary);

	cout << "There are " << voroVertices.size() << " Voronoi vertices" << endl;
	cout << "There are " << delaunayVertices.size() << " Delaunay vertices" << endl;
	vector<Vector3D> diff;
	set_difference(voroVertices.begin(), voroVertices.end(), delaunayVertices.begin(), delaunayVertices.end(), back_inserter(diff));

	cout << endl << "The difference is: " << endl;
	for (auto v : diff)
	{
		bool inside = boundary.inside(v);
		cout << namer.GetName(v) << v << ": " << (inside ? "inside" : "outside") << endl;
	}

/*	cout << "\nVertices:\n==========\n";
	for(auto pair: namer)
		cout << pair.first << ": " << pair.second << endl; */

#ifdef _DEBUG
	string s;
	getline(cin, s);
#endif

	return 0;
}

set<Vector3D> UseVoroPlusPlus(const vector<Vector3D>& points, const OuterBoundary3D &boundary)
{
	VoroPlusPlus tes;
	tes.Initialise(points, boundary);

	cout << "Voronoi:\n========" << endl;
	set<Vector3D> vertices;
	for (unsigned i = 0; i < tes.GetPointNo(); i++)
	{
		cout << "Cell " << namer.GetName(points[i]) << " at " << points[i] << endl;
		auto faces = tes.GetCellFaces(i);
		for (unsigned j = 0; j < faces.size(); j++)
		{
			auto face = tes.GetFace(faces[j]);
			cout << "\tFace " << j << " with neighbor " << face.OtherNeighbor(i) << endl << "\t\t";
			for (auto it = face.vertices.begin(); it != face.vertices.end(); it++)
			{
				cout << namer.GetName(*it) << " ";
				vertices.insert(*it);
			}
			cout << endl;
		}
	}

	return vertices;
}

Tetrahedron FindBigTetrahedron(const OuterBoundary3D &boundary)
{
	// A big tetrahedron that will contain the bounding box, as well as the 8 adjacent boundary boxes,
	// and with room to spare.

	const Vector3D fur = boundary.FrontUpperRight();
	const Vector3D bll = boundary.BackLowerLeft();
	Vector3D absFrontUpperRight(abs(fur.x), abs(fur.y), abs(fur.z));
	Vector3D absBackLowerLeft(abs(bll.x), abs(bll.y), abs(bll.z));

	absFrontUpperRight *= 1000;
	absBackLowerLeft *= -1000;

	// The top of the tetrahedron will be on the Y axis
	vector<Vector3D> tetrahedron;
	tetrahedron.push_back(Vector3D(0, absFrontUpperRight.y, 0));

	// The bottom face is parallel to the x-z plane
	double bottomY = absBackLowerLeft.y;

	// The bottom face is a triangle whose lower edge is parallel to the x axis
	double backZ = absBackLowerLeft.z;
	tetrahedron.push_back(Vector3D(absBackLowerLeft.x, bottomY, backZ));
	tetrahedron.push_back(Vector3D(absFrontUpperRight.x, bottomY, backZ));

	// The last triangle edge is on x=0
	tetrahedron.push_back(Vector3D(0, bottomY, absFrontUpperRight.z));

	return Tetrahedron(tetrahedron);
}

//\brief Find all the Outer Tetrahedra
//\remark An Outer Tetrahedron is a tetrahedron that has a vertex in the big tetrahedron
unordered_set<int> FindOuterTetrahedra(const Delaunay &del)
{
	unordered_set<int> result;

	for (int i = 0; i < 4; i++)
	{
		Vector3D pt = del.BigTetrahedron()[i];
		vector<int> tetrahedra = del.VertexNeighbors(pt);
		result.insert(tetrahedra.begin(), tetrahedra.end());
	}

	return result;
}

//\brief Find all the Edge Tetrahedra
//\remark An Edge Tetrahedron has a vertex that belongs to an Outer Tetrahedron. Edge Tetrahedra are not Outer Tetrahedra
unordered_set<int> FindEdgeTetrahedra(const Delaunay &del, const unordered_set<int>& outerTetrahedra)
{
	unordered_set<int> edgeTetrahedra;

	// Go over all the outer tetrahedra
	for (unordered_set<int>::const_iterator it = outerTetrahedra.begin(); it != outerTetrahedra.end(); it++)
	{
		Tetrahedron t = del[*it];
		// And all their vertices
		for (int i = 0; i < 4; i++)
		{
			Vector3D pt = t[i];
			if (del.IsBigTetrahedron(pt))  // Ignore the BigTetrahedron - all tetrahedra touching this vertex are Outer and not Edge
				continue;
			// Find all the neighbor tetrahedra of pt
			vector<int> tetrahedra = del.VertexNeighbors(pt);
			// Add them to the set
			for (vector<int>::iterator ptIt = tetrahedra.begin(); ptIt != tetrahedra.end(); ptIt++)
				if (!contains(outerTetrahedra, *ptIt))  // But only if their not Outer
					edgeTetrahedra.insert(*ptIt);
		}
	}

	return edgeTetrahedra;
}


set<Vector3D> UseTetGen(const vector<Vector3D>& points, const OuterBoundary3D &boundary)
{
	cout << "TetGen:\n=======\n";

	auto bigTetrahedron = FindBigTetrahedron(boundary);
	cout << "Big Tetrahedron is " << bigTetrahedron << endl;
	for (int i = 0; i < 4; i++)
		namer.GetName(bigTetrahedron[i], "B");

	cout << "Running Delaunay for the first time" << endl;
	TetGenDelaunay tetgen1(points, bigTetrahedron);
	tetgen1.Run();
	cout << "Delaunay returned " << tetgen1.NumTetrahedra() << " tetrahedra" << endl;

	cout << "Calculating ghost points..." << endl;
	//auto ghosts = GetGhostPoints(tetgen1, boundary);
	auto ghosts = BruteForceGhosts(tetgen1, boundary);

	cout << "There are " << ghosts.size() << " ghost points" << endl;
	cout << "Running the second Delaunay..." << endl;
	vector<Vector3D> allPoints(points);
	allPoints.insert(allPoints.end(), ghosts.begin(), ghosts.end());
	TetGenDelaunay tetgen2(allPoints, bigTetrahedron);
	tetgen2.Run();
	cout << "Second Delaunay returned " << tetgen2.NumTetrahedra() << " tetrahedra" << endl;

	set<Vector3D> realVertices;
	for (int i = 0; i < tetgen2.NumTetrahedra(); i++)
	{
		Tetrahedron t = tetgen2[i];
		Vector3D center = t.center();
		string name = namer.GetName(center, "D");
		cout << "Tetrahedron " << i << " centered at " << name << center << endl;
		if (boundary.inside(center))
		{
			realVertices.insert(center);
		}
	}

	return realVertices;
}

map<Vector3D, unordered_set<Subcube>> FindHullBreaches(const Delaunay &del, const unordered_set<int>& edgeTetrahedra, 
	const unordered_set<int> &outerTetrahedra, const OuterBoundary3D &boundary)
{
	map<Vector3D, unordered_set<Subcube>> result;
	for (unordered_set<int>::iterator itTetra = edgeTetrahedra.begin(); itTetra != edgeTetrahedra.end(); itTetra++)
	{
		Tetrahedron t = del[*itTetra];
		for (int iv = 0; iv < 4; iv++)
		{
			Vector3D pt = t[iv];
			if (result.find(pt) != result.end())
				continue;

			vector<int> tetrahedraIndices = del.VertexNeighbors(pt);
			unordered_set<Subcube> breaches;

			for (vector<int>::iterator itTetrahedron = tetrahedraIndices.begin(); itTetrahedron != tetrahedraIndices.end(); itTetrahedron++)
			{
				if (contains(outerTetrahedra, *itTetrahedron))
					continue;

				Tetrahedron tetrahedron = del[*itTetrahedron];
				for (set<Subcube>::iterator itSubcube = Subcube::all().begin(); itSubcube != Subcube::all().end(); itSubcube++)
				{
					if (boundary.distance(tetrahedron.center(), *itSubcube) < tetrahedron.radius())
						breaches.insert(*itSubcube);
				}
			}
			result[pt] = breaches;
		}
	}

	return result;
}

set<Vector3D> GetGhostPoints(const Delaunay &del, const OuterBoundary3D &boundary)
{
	unordered_set<int> outer = FindOuterTetrahedra(del);
	unordered_set<int> edge = FindEdgeTetrahedra(del, outer);
	map<Vector3D, unordered_set<Subcube>> breaches = FindHullBreaches(del, edge, outer, boundary);

	set<Vector3D> ghosts;
	if (breaches.empty())
		return ghosts;

	for (map<Vector3D, unordered_set<Subcube>>::iterator it = breaches.begin(); it != breaches.end(); it++)
	{
		Vector3D pt = it->first;
		for (unordered_set<Subcube>::iterator itSubcube = it->second.begin(); itSubcube != it->second.end(); itSubcube++)
		{
			Vector3D ghost = boundary.ghost(pt, *itSubcube);
			ghosts.insert(ghost);
		}
	}

	return ghosts;
}

set<Vector3D> BruteForceGhosts(const Delaunay &del, const OuterBoundary3D &boundary)
{
	set<Vector3D> ghosts;
	for (auto vec : del.InputPoints())
		for (auto subcube : Subcube::all())
			ghosts.insert(boundary.ghost(vec, subcube));

	return ghosts;
}