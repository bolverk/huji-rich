// StressVoronoi.cpp : Defines the entry point for the console application.
//

#include "GeometryCommon/Vector3D.hpp"
#include "GeometryCommon/OuterBoundary3D.hpp"
#include "Voronoi/VoroPlusPlus.hpp"
#include "Voronoi/DelaunayVoronoi.hpp"
#include "Voronoi/GhostBusters.hpp"
#include "Voronoi/TetGenDelaunay.hpp"

#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

using namespace std;

void InitRandom();
double RandomDouble(double min, double max);
Vector3D RandomPoint(const OuterBoundary3D &boundary);
void RandomBoundary(Vector3D &frontTopRight, Vector3D &backLowerLeft);
vector<Vector3D> RandomPoints(int num, const OuterBoundary3D &boundary);

void WritePoints(const vector<Vector3D> &points, const OuterBoundary3D &boundary, const std::string &filename);

void RunVoronoi(const vector<Vector3D> &points, const OuterBoundary3D &boundary, Tessellation3D &tes);
bool CompareTessallations(const Tessellation3D &tes1, const Tessellation3D &tes2);

void TestDelaunay(const vector<Vector3D> &points, const OuterBoundary3D &boundary);

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

	string GetName(VectorRef vec, string prefix = "V")
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

VectorNamer namer;

int main(int argc, char*argv[])
{
	namer.GetName(Vector3D(0, 0, 0), "Z");

	InitRandom();

	Vector3D ftr, bll;
	RandomBoundary(ftr, bll);
	RectangularBoundary3D boundary(ftr, bll);

	vector<Vector3D> points = RandomPoints(7, boundary);
	WritePoints(points, boundary, "stress.node");

	cout << "Generated " << points.size() << " points" << endl;

	cout << "Testing Delaunay...";
	TestDelaunay(points, boundary);

	cout << "Running voro++..." << endl;
	VoroPlusPlus voro;
	RunVoronoi(points, boundary, voro);

	cout << "Running tetgen..." << endl;
	DelaunayVoronoi<TetGenDelaunay, BruteForceGhostBuster>  del;
	RunVoronoi(points, boundary, del);

	cout << "Comparing results..." << endl;
	bool same = CompareTessallations(voro, del);
	if (same)
		cout << "Success!!!" << endl;
	else
		cout << "Difference" << endl;

#ifdef _DEBUG
	string s;
	getline(cin, s);
#endif
}

void InitRandom()
{
	srand(17);
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
	backLowerLeft.x = -1000;
	frontTopRight.x = 1000;

	backLowerLeft.y = -1000;
	frontTopRight.y = 1000;
	
	backLowerLeft.z = -1000;
	frontTopRight.z = 1000;
}

vector<Vector3D> RandomPoints(int num, const OuterBoundary3D &boundary)
{
	vector<Vector3D> points(num);

	for (int i = 0; i < num; i++)
		points[i] = RandomPoint(boundary);

	return points;
}

void RunVoronoi(const vector<Vector3D> &points, const OuterBoundary3D &boundary, Tessellation3D &tes)
{
	tes.Initialise(points, boundary);
}

bool CompareTessallations(const Tessellation3D &tes1, const Tessellation3D &tes2)
{
	if (tes1.GetPointNo() != tes2.GetPointNo())
	{
		cout << "Not the same point number" << endl;
		return false;
	}

	for (int i = 0; i < tes1.GetPointNo(); i++)
	{
		double volume1 = tes1.GetVolume(i);
		double volume2 = tes2.GetVolume(i);
		if (volume1 != volume2)
			cout << "Different volume in cell " << i << endl;

		vector<size_t> faces1 = tes1.GetCellFaces(i);
		vector<size_t> faces2 = tes2.GetCellFaces(i);
		if (faces1.size() != faces2.size())
			cout << "Wrong number of faces in cell " << i << endl;
	}

	return true;
}

void WritePoints(const vector<Vector3D> &points, const OuterBoundary3D &boundary, const std::string &filename)
{
	Tetrahedron big = Delaunay::CalcBigTetrahedron(boundary);
	ofstream output;
	
	output.open(filename);
	output << "# Nodes, Dim, #Attrs, Boundary" << endl;
	output << points.size() + 4 << ",3,0,1" << endl;
	output << "# Index, X, Y, Z, Boundary" << endl;
	output << "# Main points" << endl;

	for (int i = 0; i < points.size(); i++)
		output << i + 1 << " " << points[i].x << " " << points[i].y << " " << points[i].z << " 0" << endl;

	output << "# Boundary points" << endl;
	for (int i = 0; i < 4; i++)
		output << i + 1 + points.size() << " " << big[i]->x << " " << big[i]->y << " " << big[i]->z << " 1" << endl;

	output << endl;
}


void TestDelaunay(const vector<Vector3D> &points, const OuterBoundary3D &boundary)
{
	for (auto point : points)
		namer.GetName(point, "C");
	Tetrahedron big = Delaunay::CalcBigTetrahedron(boundary);

	for (int i = 0; i < 4; i++)
		namer.GetName(big[i], "B");

	TetGenDelaunay tetgen(VectorRef::vector(points), big);
	tetgen.Run();

	cout << "There are " << tetgen.NumTetrahedra() << " resulting tetrahedra" << endl;
	for (int i = 0; i < tetgen.NumTetrahedra(); i++)
	{
		VectorRef center = tetgen[i].center();
		cout << "\t" << i + 1 << " " << namer.GetName(center) << center << ": ";
		cout << "\t\t";
		for (int j = 0; j < 4; j++)
			cout << namer.GetName(tetgen[i][j]) << " ";
		cout << endl;
	}
}
