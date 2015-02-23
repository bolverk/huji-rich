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
#include <fstream>
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
void WriteTetGenInput(const vector<Vector3D> &meshPoints, const OuterBoundary3D &boundary,  const vector<VectorRef> &allPoints, const std::string &folder);
//void CalcCenters(const vector<Vector3D> &points, const OuterBoundary3D &boundary);

VectorNamer namer;

int main()
{
	InitRandom();
	RectangularBoundary3D boundary(Vector3D(200, 200, 200), Vector3D(-200, -200, -200));
	/*	vector<Vector3D> vertices{ Vector3D(90, 34, 89),
		Vector3D(21, 3, 78),
		Vector3D(76, 35, 74),
		Vector3D(28, 4, 7),
		Vector3D(65, 60, 22),
		Vector3D(59, 92, 5),
		Vector3D(84, 0, 32),
		Vector3D(12, 15, 37)} */
	vector<Vector3D> vertices{
		Vector3D(-195.544, -185.156, -133.897),
		Vector3D(3.91247, 99.1302, -62.6179),
		Vector3D(101.498, -179.223, 83.3216),
		Vector3D(-131.285, 90.9024, 196.802),
		Vector3D(38.9477, 1.58086, -59.2608),
		Vector3D(-40.9864, -154.003, -187.414),
		Vector3D(71.3584, -104.16, -9.00296),
		Vector3D(-113.071, -7.00095, -121.018),
		Vector3D(-148.79, 115.647, 163.976),
		Vector3D(-60.7257, 117.161, -121.86),
	};

	DelaunayVoronoi<TetGenDelaunay, BruteForceGhostBuster> del;
	del.Initialise(vertices, boundary);
	WriteTetGenInput(vertices, boundary, del.AllPoints, "..\\PythonVoronoi\\points\\");

	/*auto vertices = RandomPoints(10, boundary);
	for (auto vertice : vertices)
		cout << vertice << endl; */

	namer.GetName(Vector3D(), "Z"); // Zero vector
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
		cout << "Cell " << namer.GetName(points[i]) << " at " << points[i] << " with " << faces.size() << " faces" << endl;
		// cout << "\tVolume: " << tes.GetVolume(i) << ", Center of Mass: " << namer.GetName(tes.GetCellCM(i)) << tes.GetCellCM(i) << endl;
		for (unsigned j = 0; j < faces.size(); j++)
		{
			auto face = tes.GetFace(faces[j]);

			cout << "\tFace F" << faces[j];
			if (face.OtherNeighbor(i).is_initialized())
				cout << " neighbor C" << face.OtherNeighbor(i).value().GetCell() + 1;
			cout << endl;
			for (auto it = face.vertices.begin(); it != face.vertices.end(); it++)
			{
				cout << "\t\t" << namer.GetName(*it) << " " << *it << endl;
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

void WriteTetGenInput(const vector<Vector3D> &points, const OuterBoundary3D &boundary, const vector<VectorRef>& allPoints, const std::string &folder)
{
	// First, write the original points
	ofstream origOutput;
	origOutput.open(folder + "orig.node");

	origOutput << points.size() + 4 << " 3 0 1" << endl;
	origOutput << "# Mesh points" << endl;
	int index = 1;
	for (auto point : points)
		origOutput << index++ << " " << point.x << " " << point.y << " " << point.z << " 0" << endl;

	auto big = Delaunay::CalcBigTetrahedron(boundary);
	for (int i = 0; i < 4; i++)
		origOutput << index++ << " " << big[i]->x << " " << big[i]->y << " " << big[i]->z << " 1" << endl;

	origOutput << endl << "# Generated automatically by TryVoronois" << endl;

	// Now write *all* the points, including ghosts
	ofstream allOutput;
	allOutput.open(folder + "all.node");
	allOutput << allPoints.size() + 4 << " 3 0 1" << endl;
	index = 1;
	for (auto point: allPoints)
		allOutput << index++ << " " << point->x << " " << point->y << " " << point->z << " 0" << endl;
	for (int i = 0; i < 4; i++)
		allOutput << index++ << " " << big[i]->x << " " << big[i]->y << " " << big[i]->z << " 1" << endl;
	allOutput << endl << "# Generated automatically by TryVoronois" << endl;
}

void CalcCenters(const vector<Vector3D> &points, const OuterBoundary3D &boundary)
{
	auto big = Delaunay::CalcBigTetrahedron(boundary);
	TetGenDelaunay tet(VectorRef::vector(points), big);
	tet.Run();

	for (size_t i = 0; i < tet.NumTetrahedra(); i++)
	{
		cout << i << ": " << *tet[i].center() << endl;
	}
}
