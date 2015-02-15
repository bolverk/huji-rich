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

void InitRandom();
double RandomDouble(double min, double max);
Vector3D RandomPoint(const OuterBoundary3D &boundary);
void RandomBoundary(Vector3D &frontTopRight, Vector3D &backLowerLeft);

vector<Vector3D> RandomPoints(int num, const OuterBoundary3D &boundary);

class VectorNamer
{
private:
	map<VectorRef, string> _vectors;
	int _last;

public:
	VectorNamer()
	{
		_last = 0;
	}

	string GetName(VectorRef vec, string prefix="V")
	{
		auto inMap = _vectors.find(vec);
		if (inMap != _vectors.end())
			return inMap->second;

		stringstream strm;
		strm << prefix << _last++;
		_vectors[vec] = strm.str();
		return strm.str();
	}

	map<VectorRef, string>::const_iterator begin() const { return _vectors.cbegin(); }
	map<VectorRef, string>::const_iterator end() const { return _vectors.cend(); }
};

void UseVoroPlusPlus(const vector<Vector3D>& points, const OuterBoundary3D &boundary);
void UseTetGen(const vector<Vector3D>& points, const OuterBoundary3D &boundary);

void DisplayResults(const vector<Vector3D> &points, const Tessellation3D &tes);

VectorNamer namer;

int main()
{
	InitRandom();
	RectangularBoundary3D boundary(Vector3D(200, 200, 200), Vector3D(-200, -200, -200));
	vector<Vector3D> vertices{ Vector3D(90, 34, 89),
		Vector3D(21, 3, 78),
		Vector3D(76, 35, 74),
		Vector3D(28, 4, 7),
		Vector3D(65, 60, 22),
		Vector3D(59, 92, 5),
		Vector3D(84, 0, 32),
		Vector3D(12, 15, 37),
		Vector3D(-17, -17, -17)};
	//auto vertices = RandomPoints(7, boundary);

	for (auto vec : vertices)
		namer.GetName(vec, "C");

	UseVoroPlusPlus(vertices, boundary);
	UseTetGen(vertices, boundary);

#ifdef _DEBUG
	string s;
	getline(cin, s);
#endif

	return 0;
}

void UseVoroPlusPlus(const vector<Vector3D>& points, const OuterBoundary3D &boundary)
{
	VoroPlusPlus tes;
	tes.Initialise(points, boundary);

	cout << "Voronoi:\n========" << endl;
	DisplayResults(points, tes);
}

void DisplayResults(const vector<Vector3D> &points, const Tessellation3D &tes)
{
	set<Vector3D> vertices;
	for (unsigned i = 0; i < tes.GetPointNo(); i++)
	{
		auto faces = tes.GetCellFaces(i);
		cout << "Cell " << namer.GetName(points[i]) << " at " << points[i] << " with " << faces.size() << "faces" << endl;
		cout << "\tVolume: " << tes.GetVolume(i) << ", Center of Mass: " << namer.GetName(tes.GetCellCM(i)) << tes.GetCellCM(i) << endl;
		for (unsigned j = 0; j < faces.size(); j++)
		{
			auto face = tes.GetFace(faces[j]);
			cout << "\tFace F" << faces[j] << " with neighbor cell " << face.OtherNeighbor(i) << endl << "\t\t";
			for (auto it = face.vertices.begin(); it != face.vertices.end(); it++)
			{
				cout << namer.GetName(*it) << " ";
				vertices.insert(**it);
			}
			cout << endl;
		}
	}
}

void UseTetGen(const vector<Vector3D>& points, const OuterBoundary3D &boundary)
{
	cout << "TetGen:\n=======\n";
	DelaunayVoronoi<TetGenDelaunay, BruteForceGhostBuster> delaunay;
	delaunay.Initialise(points, boundary);

	DisplayResults(points, delaunay);
}

void InitRandom()
{
	srand(100);
}

double RandomDouble(double min, double max)
{
	double fraction = (double)rand() / RAND_MAX;
	return (max - min) * fraction + min;
}

Vector3D RandomPoint(const OuterBoundary3D &boundary)
{
	Vector3D v;

	v.x = RandomDouble(boundary.BackLowerLeft().x, boundary.FrontUpperRight().x);
	v.y = RandomDouble(boundary.BackLowerLeft().y, boundary.FrontUpperRight().y);
	v.z = RandomDouble(boundary.BackLowerLeft().z, boundary.FrontUpperRight().z);

	return v;
}

void RandomBoundary(Vector3D &frontTopRight, Vector3D &backLowerLeft)
{
	backLowerLeft.x = 100;
	frontTopRight.x = 1000;

	backLowerLeft.y = -500;
	frontTopRight.y = 500;

	backLowerLeft.z = -1500;
	frontTopRight.z = -100;
}

vector<Vector3D> RandomPoints(int num, const OuterBoundary3D &boundary)
{
	vector<Vector3D> points(num);

	for (int i = 0; i < num; i++)
		points[i] = RandomPoint(boundary);

	return points;
}

