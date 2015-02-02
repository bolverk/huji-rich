// Voronoi.cpp : Defines the entry point for the console application.
//

#include "tetgen.h"
#include "Mat44.hpp"
#include "GeometryCommon/Vector3D.hpp"
#include "Voronoi/VoroPlusPlus.hpp"
#include <vector>
#include <iostream>
#include <set>
#include "Tetrahedron.hpp"

using namespace std;

void UseVoroPlusPlus(const vector<Vector3D>& points);

void RunTetGen(const vector<Vector3D> &points, tetgenio &out);
void UseTetGen(const vector<Vector3D>& points);
Tetrahedron FindBigTetrahedron(const OuterBoundary3D &boundary);

int main()
{
	vector<Vector3D> vertices{ Vector3D(90, 34, 89),
		Vector3D(21, 3, 78),
		Vector3D(76, 35, 74),
		Vector3D(28, 4, 7),
		Vector3D(65, 60, 22),
		Vector3D(59, 92, 5),
		Vector3D(84, 0, 32) };

	UseVoroPlusPlus(vertices);
	UseTetGen(vertices);

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

	for (unsigned i = 0; i < tes.GetPointNo(); i++)
	{
		auto faces = tes.GetCellFaces(i);
		for (unsigned j = 0; j < faces.size(); j++)
		{
			auto face = tes.GetFace(j);
			for (auto it = face.vertices.begin(); it != face.vertices.end(); it++)
			{
				//if (it->x > -200 && it->x<200 && it->y > -200 && it->y < 200 && it->z > -200 && it->z < 200)
					vertices.insert(*it);
			}
		}
	}

	cout << "There are " << vertices.size() << " non-infinity vertices" << endl;
	for (auto it = vertices.begin(); it != vertices.end(); it++)
		cout << *it << endl;
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
	b.parse_commandline("efnQ");
	tetrahedralize(&b, &in, &out);
}

void UseTetGen(const vector<Vector3D>& points)
{
	OuterBoundary3D boundary(OuterBoundary3D::RECTANGULAR, Vector3D(200, 200, 200), Vector3D(-200, -200, -200));
	tetgenio out;

	RunTetGen(points, out);
	cout << "Calculated " << out.numberoftetrahedra << " tetrahedra" << endl;
	vector<Tetrahedron> tetrahedra;
	int offset = 0;
	for (int i = 0; i < out.numberoftetrahedra; i++)
	{
		vector<Vector3D> vertices;
		for (int j = 0; j < 4; j++)
		{
			int ptIndex = out.tetrahedronlist[offset++];
			vertices.push_back(points[ptIndex]);
		}
		Tetrahedron t(vertices);
		tetrahedra.push_back(t);
		// cout << t << ": " << t.center() << endl;
		cout << t.center() << endl;
	}
}

