// Voronoi.cpp : Defines the entry point for the console application.
//

#include "tetgen.h"
#include "Mat44.hpp"
#include "GeometryCommon/Vector3D.hpp"
#include "Voronoi/VoroPlusPlus.hpp"
#include <vector>
#include <iostream>
#include <set>

using namespace std;

void UseVoroPlusPlus(const vector<Vector3D>& points);
void UseTetGen(const vector<Vector3D>& points);
Vector3D FindCircumcenter(const vector<Vector3D>& vertices);

int main()
{
	vector<Vector3D> vertices{ Vector3D(90, 34, 89),
		Vector3D(21, 3, 78),
		Vector3D(76, 35, 74),
		Vector3D(28, 4, 7),
		Vector3D(65, 60, 22),
		Vector3D(59, 92, 5),
		Vector3D(84, 0, 32) };

	// UseVoroPlusPlus(vertices);
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

void UseTetGen(const vector<Vector3D>& points)
{
	tetgenio in, out;

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
	b.parse_commandline("-efn");
	tetrahedralize(&b, &in, &out);

	cout << "Calculated " << out.numberoftetrahedra << " tetrahedra" << endl;
	vector<vector<int>> tetrahedra;
	int offset = 0;
	for (int i = 0; i < out.numberoftetrahedra; i++)
	{
		vector<int> indices(4);
		for (int j = 0; j < 4; j++)
			indices[j] = out.tetrahedronlist[offset++];
		tetrahedra.push_back(indices);
	}

	for (int i = 0; i < tetrahedra.size(); i++)
	{
		vector<Vector3D> vertices(4);
		for (int j = 0; j < 4; j++)
		{
			vertices[j] = points[tetrahedra[i][j]];
			cout << vertices[j] << " ";
		}
		Vector3D center = FindCircumcenter(vertices);
		cout << ": " << center << endl;
	}
}

// \brief Find the circumcenter of a tetrahedron
// \param vertices - a vector of the 4 corners
// \returns The circumcenter
// \remark Taken from here: http://mathworld.wolfram.com/Circumsphere.html
Vector3D FindCircumcenter(const vector<Vector3D>& vertices)
{
	if (vertices.size() != 4)
		throw invalid_argument("Only tetrahedra are supported");

	Vector3D v1 = vertices[0];
	Vector3D v2 = vertices[1];
	Vector3D v3 = vertices[2];
	Vector3D v4 = vertices[3];

	
	Mat44<double> m_a{ v1.x, v1.y, v1.z, 1,
		v2.x, v2.y, v2.z, 1,
		v3.x, v3.y, v3.z, 1,
		v4.x, v4.y, v4.z, 1 };
	double a = m_a.determinant();

	Mat44<double> m_Dx = { abs2(v1), v1.y, v1.z, 1,
		abs2(v2), v2.y, v2.z, 1,
		abs2(v3), v3.y, v3.z, 1,
		abs2(v4), v4.y, v4.z, 1 };
	double Dx = m_Dx.determinant();

	Mat44<double> m_Dy = { abs2(v1), v1.x, v1.z, 1,
		abs2(v2), v2.x, v2.z, 1,
		abs2(v3), v3.x, v3.z, 1,
		abs2(v4), v4.x, v4.z, 1 };
	double Dy = -m_Dy.determinant();

	Mat44<double> m_Dz = { abs2(v1), v1.x, v1.y, 1,
		abs2(v2), v2.x, v2.y, 1,
		abs2(v3), v3.x, v3.y, 1,
		abs2(v4), v4.x, v4.y, 1 };
	double Dz = m_Dz.determinant();

	Mat44<double> m_c = { abs2(v1), v1.x, v1.y, v1.z,
		abs2(v2), v2.x, v2.y, v2.z,
		abs2(v3), v3.x, v3.y, v3.z,
		abs2(v4), v4.x, v4.y, v4.z };
	double c = m_c.determinant();

	return Vector3D(Dx / (2 * a), Dy / (2 * a), Dz / (2 * a));
}

