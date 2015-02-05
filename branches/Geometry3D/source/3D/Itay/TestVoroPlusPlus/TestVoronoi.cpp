// TestVoroPlusPlus.cpp : Defines the entry point for the console application.
//

#include "gtest/gtest.h"
#include "GeometryCommon/Face.hpp"
#include "GeometryCommon/Vector3D.hpp"
#include "Voronoi/TessellationBase.hpp"
#include "Voronoi/VoroPlusPlus.hpp"
#include "Utilities/assert.hpp"
#include "GeometryCommon/OuterBoundary3D.hpp"

#include <vector>
#include <cstdlib>
#include <cassert>

#ifdef _DEBUG
#include <iostream>
#include <string>
#endif

using namespace std;

TEST(VoroPlusPlus, FaceStore)
{
	TessellationBase::FaceStore store;
	vector<Face> faces;

	for (int i = 0; i < 100; i++)
	{
		vector<Vector3D> vertices;
		int numVertices = rand() % 6 + 2;
		for (int j = 0; j < numVertices; j++)
			vertices.push_back(Vector3D(rand(), rand(), rand()));
		size_t index = store.StoreFace(vertices);
		ASSERT_EQ(index, i);
		faces.push_back(store.GetFace(index));
	}

	for (int i = 99; i >= 0; i--)
	{
		size_t index = store.StoreFace(faces[i].vertices);
		ASSERT_EQ(index, i);
	}
}


double RandomDouble(double min, double max)
{
	double f = (double)rand() / RAND_MAX;  // Between 0 and 1,inclusive
	return (max - min) * f + min;
}

Vector3D RandomVector(const OuterBoundary3D &boundary)
{
	Vector3D v;

	v.x = RandomDouble(boundary.BackLowerLeft().x, boundary.FrontUpperRight().x);
	v.y = RandomDouble(boundary.BackLowerLeft().y, boundary.FrontUpperRight().y);
	v.z = RandomDouble(boundary.BackLowerLeft().z, boundary.FrontUpperRight().z);

	return v;
}

vector<Vector3D> CreateRandom(int numPoints, const OuterBoundary3D &boundary)
{
	srand(17);   // Make sure we generate the same sequence over and over
	vector<Vector3D> points;

	for (int i = 0; i < numPoints; i++)
	{
		Vector3D v = RandomVector(boundary);
		points.push_back(v);
	}

	return points;
}

vector<Vector3D> CreateCubeMesh(int perSide)
{
	vector<Vector3D> mesh;
	for (int x = 0; x < perSide; x++)
		for (int y = 0; y < perSide; y++)
			for (int z = 0; z < perSide; z++)
				mesh.push_back(Vector3D(x, y, z));

	return mesh;
}

void EnsureCubeVoronoi(const vector<Vector3D> &mesh, const Tessellation3D &tes)
{
	ASSERT_EQ(tes.GetPointNo(), mesh.size());
	for (int pt = 0; pt < mesh.size(); pt++)
	{
		auto faces = tes.GetCellFaces(pt);
		EXPECT_EQ(faces.size(), 6);
		for (int fc = 0; fc < faces.size(); fc++)
		{
			auto face = tes.GetFace(fc);
			EXPECT_EQ(face.vertices.size(), 4);
			EXPECT_NEAR(face.GetArea(), 1.0, 1e-12);
		}
		auto CoM = tes.GetCellCM(pt);
		EXPECT_EQ(CoM, mesh[pt]);
		EXPECT_NEAR(tes.GetVolume(pt), 1.0, 1e-12);
	}
}

TEST(VoroPlusPlus, Cube)
{
	vector<Vector3D> mesh = CreateCubeMesh(5);

	VoroPlusPlus tes;
	RectangularBoundary3D boundary(Vector3D(4.5, 4.5, 4.5),
		Vector3D(-0.5, -0.5, -0.5));
	tes.Initialise(mesh, boundary);
	EnsureCubeVoronoi(mesh, tes);
}

void assertion_gtest_bridge(const char *expr, const char *function, const char *file, long line)
{
	EXPECT_FALSE(expr);
}

int main(int argc, char *argv[])
{
	BOOST_ASSERT_HANDLER = assertion_gtest_bridge;
	testing::InitGoogleTest(&argc, argv);
	int rc = RUN_ALL_TESTS();

#ifdef _DEBUG
	string s;
	getline(cin, s);
#endif
}