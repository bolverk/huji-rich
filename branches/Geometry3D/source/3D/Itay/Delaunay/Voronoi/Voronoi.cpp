// Voronoi.cpp : Defines the entry point for the console application.
//

#include "tetgen.h"
#include "GeometryCommon/Mat44.hpp"
#include "GeometryCommon/Vector3D.hpp"
#include "Voronoi/VoroPlusPlus.hpp"
#include <vector>
#include <iostream>
#include <set>
#include "GeometryCommon/Tetrahedron.hpp"
#include <map>
#include <sstream>
#include <string>

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

void UseVoroPlusPlus(const vector<Vector3D>& points);

void RunTetGen(const vector<Vector3D> &points, tetgenio &out);
void UseTetGen(const vector<Vector3D>& points);
Tetrahedron FindBigTetrahedron(const OuterBoundary3D &boundary);

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
	for (auto vec : vertices)
		namer.GetName(vec, "C");

	UseVoroPlusPlus(vertices);
	UseTetGen(vertices);

	cout << "\nVertices:\n==========\n";
	for(auto pair: namer)
		cout << pair.first << ": " << pair.second << endl;

#ifdef _DEBUG
	string s;
	getline(cin, s);
#endif

	return 0;
}

void UseVoroPlusPlus(const vector<Vector3D>& points)
{
	VoroPlusPlus tes;
	OuterBoundary3D boundary(OuterBoundary3D::RECTANGULAR, Vector3D(200, 200, 200), Vector3D(-200, -200, -200));
	tes.Initialise(points, boundary);
	set<Vector3D> vertices;

	cout << "Voronoi:\n========" << endl;
	for (unsigned i = 0; i < tes.GetPointNo(); i++)
	{
		cout << "Cell " << namer.GetName(points[i]) << " at " << points[i] << endl;
		auto faces = tes.GetCellFaces(i);
		for (unsigned j = 0; j < faces.size(); j++)
		{
			auto face = tes.GetFace(faces[j]);
			cout << "\tFace " << j << " with neighbor " << face.OtherNeighbor(i) << endl << "\t\t";
			for (auto it = face.vertices.begin(); it != face.vertices.end(); it++)
				cout << namer.GetName(*it) << " ";
			cout << endl;
		}
	}
}

Tetrahedron FindBigTetrahedron(const OuterBoundary3D &boundary)
{
	// A big tetrahedron that will contain the bounding box, as well as the 8 adjacent boundary boxes,
	// and with room to spare.

	const Vector3D fur = boundary.FrontUpperRight();
	const Vector3D bll = boundary.BackLowerLeft();
	Vector3D absFrontUpperRight(abs(fur.x), abs(fur.y), abs(fur.z));
	Vector3D absBackLowerLeft(abs(bll.x), abs(bll.y), abs(bll.z));

	absFrontUpperRight *= 10;
	absBackLowerLeft *= -10;

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

void RunTetGen(const vector<Vector3D> &points, tetgenio &out)
{
	tetgenio in;

	in.firstnumber = 0;
	in.pointlist = new REAL[points.size() * 3];
	int index = 0;
	for (auto it = points.begin(); it != points.end(); it++)
	{
		in.pointlist[index++] = it->x;
		in.pointlist[index++] = it->y;
		in.pointlist[index++] = it->z;
	}
	in.numberofpoints = (int)points.size();

	tetgenbehavior b;
	b.parse_commandline("efnvQ");
	tetrahedralize(&b, &in, &out);
}

void UseTetGen(const vector<Vector3D>& points)
{
	OuterBoundary3D boundary(OuterBoundary3D::RECTANGULAR, Vector3D(200, 200, 200), Vector3D(-200, -200, -200));
	tetgenio out;

	cout << "TetGen:\n=======" << endl;

	RunTetGen(points, out);
	vector<Tetrahedron> tetrahedra;
	vector<vector<int>> tetrahedra_neighbors;

	out.vedgelist;
	
	int offset = 0;
	for (int i = 0; i < out.numberoftetrahedra; i++)
	{
		vector<Vector3D> vertices;
		vector<int> neighbors;
		for (int j = 0; j < 4; j++)
		{
			int ptIndex = out.tetrahedronlist[offset];
			vertices.push_back(points[ptIndex]);

			neighbors.push_back(out.neighborlist[offset]);
			offset++;
		}
		Tetrahedron t(vertices);
		tetrahedra.push_back(t);
		tetrahedra_neighbors.push_back(neighbors);

		cout << "Vertex " <<i << ": " << namer.GetName(t.center()) << "(" << t.center() << ")" << endl << "\t";
		for (auto v : t.vertices())
			cout << namer.GetName(v) << " ";
		cout << endl << "\t";
		for (int n : neighbors)
			cout << n << " ";
		cout << endl;
	}
}

