// TestVoroPlusPlus.cpp : Defines the entry point for the console application.
//

#include "gtest/gtest.h"
#include "GeometryCommon/Face.hpp"
#include "GeometryCommon/Vector3D.hpp"
#include "Voronoi/VoroPlusPlus.hpp"
#include "Utilities/assert.hpp"

#include <vector>
#include <cstdlib>
#include <cassert>

#ifdef _DEBUG
#include <iostream>
#include <string>
#endif

using namespace std;

TEST(Geometry3D, Vector3D_Construction)
{
	Vector3D vec0;
	ASSERT_EQ(vec0.x, 0.0);
	ASSERT_EQ(vec0.y, 0.0);
	ASSERT_EQ(vec0.z, 0.0);

	Vector3D vec1(1, 2, 3);
	ASSERT_EQ(vec1.x, 1);
	ASSERT_EQ(vec1.y, 2);
	ASSERT_EQ(vec1.z, 3);

	Vector3D vec2(vec1);
	ASSERT_EQ(vec1.x, vec2.x);
	ASSERT_EQ(vec1.y, vec2.y);
	ASSERT_EQ(vec1.z, vec2.z);
}

TEST(Geometry3D, Vector3D_Operations)
{
	Vector3D v1(1, 2, 3);
	Vector3D v2(6, 5, 4);

	auto v3 = v1 + v2;
	ASSERT_EQ(v3.x, 7);
	ASSERT_EQ(v3.y, 7);
	ASSERT_EQ(v3.z, 7);

	auto v4 = v3 - v2;
	ASSERT_EQ(v4, v1);
}

TEST(Geometry3D, Face_Basic)
{
	Vector3D v1(0, 0, 0), v2(1, 0, 0), v3(1, 1, 0), v4(0, 1, 0);
	vector<Vector3D> vertices{ v1, v2, v3, v4 };
	Face face(vertices, -1, -1);

	ASSERT_EQ(vertices, face.vertices);

	Face face2(vertices, 5, 3);
	ASSERT_EQ(vertices, face2.vertices);
	ASSERT_EQ(face2.neighbors.first, 5);
	ASSERT_EQ(face2.neighbors.second, 3);
}

TEST(Geometry3D, Face_Area)
{
	Vector3D v1(0, 0, 0), v2(1, 0, 0), v3(1, 1, 0), v4(0, 1, 0);
	vector<Vector3D> vertices{ v1, v2, v3, v4 };
	Face square(vertices, -1, -1);
	ASSERT_EQ(square.GetArea(), 1);

	Vector3D v5(0, 0, 0), v6(0, 2, 0), v7(2, 0, 0);
	vector<Vector3D> vertices2{ v5, v6, v7 };
	Face triangle(vertices2, -1, -1);
	ASSERT_EQ(triangle.GetArea(), 2);
}

TEST(Geometry3D, Face_Neighbors)
{
	Vector3D v1(0, 0, 0), v2(1, 0, 0), v3(1, 1, 0), v4(0, 1, 0);
	vector<Vector3D> vertices{ v1, v2, v3, v4 };
	Face face1(vertices), face2(vertices, 1, 2), face3(vertices, 1);

	ASSERT_TRUE(face1.AddNeighbor(1));
	ASSERT_TRUE(face1.AddNeighbor(2));
	ASSERT_TRUE(face1.AddNeighbor(2));
	ASSERT_TRUE(face1.AddNeighbor(1));
	ASSERT_DEATH(face1.AddNeighbor(3), "");
	ASSERT_DEATH(face2.AddNeighbor(3), "");
	ASSERT_TRUE(face3.AddNeighbor(1));
	ASSERT_TRUE(face3.AddNeighbor(2));
	ASSERT_DEATH(face3.AddNeighbor(3), "");
}

TEST(Geometry3D, Face_Indentical)
{
	Vector3D v1(0, 0, 0), v2(1, 0, 0), v3(1, 1, 0), v4(0, 1, 0);
	vector<Vector3D> vertices1{ v1, v2, v3, v4 }, vertices2{ v2, v3, v4, v1 }, vertices3{ v4, v3, v2, v1 }, vertices4{ v1, v2, v3 };
	Face face1(vertices1, 0, 0), face2(vertices2, 0, 0), face3(vertices3, 0, 0), face4(vertices4, 0, 0);

	ASSERT_TRUE(face1.IdenticalTo(vertices2));
	ASSERT_TRUE(face1.IdenticalTo(vertices3));
	ASSERT_FALSE(face1.IdenticalTo(vertices4));
	ASSERT_TRUE(face2.IdenticalTo(vertices3));
	ASSERT_FALSE(face2.IdenticalTo(vertices4));
	ASSERT_FALSE(face3.IdenticalTo(vertices4));
}

TEST(VoroPlusPlus, FaceStore)
{
	VoroPlusPlus::FaceStore store;
	vector<Face> faces;

	for (int i = 0; i < 100; i++)
	{
		vector<Vector3D> vertices;
		int numVertices = rand() % 6 + 2;
		for (int j = 0; j < numVertices; j++)
			vertices.push_back(Vector3D(rand(), rand(), rand()));
		int index = store.StoreFace(vertices);
		ASSERT_EQ(index, i);
		faces.push_back(store.GetFace(index));
	}

	for (int i = 99; i >= 0; i--)
	{
		int index = store.StoreFace(faces[i].vertices);
		ASSERT_EQ(index, i);
	}
}

TEST(VoroPlusPlus, Cube)
{
	const int perSide = 5;

	vector<Vector3D> mesh;
	for (int x = 0; x < perSide; x++)
		for (int y = 0; y < perSide; y++)
			for (int z = 0; z < perSide; z++)
				mesh.push_back(Vector3D(x, y, z));

	VoroPlusPlus tes;
	tes.Initialise(mesh, nullptr);

	ASSERT_EQ(tes.GetPointNo(), mesh.size());
	for (int pt = 0; pt < mesh.size(); pt++)
	{
//		cout << mesh[pt].x << ", " << mesh[pt].y << ", " << mesh[pt].z << ": " << endl;
		auto faces = tes.GetCellFaces(pt);
		EXPECT_EQ(faces.size(), 6);
		for (int fc = 0; fc < faces.size(); fc++)
		{
//			cout << "\t";
			auto face = tes.GetFace(fc);
			EXPECT_EQ(face.vertices.size(), 4);
/*			for (int i = 0; i < 4; i++)
				cout << "(" << face.vertices[i].x << ", " << face.vertices[i].y << ", " << face.vertices[i].z << ") ";
			cout << endl << "\t" << "Area : " << face.GetArea() << endl; */
			EXPECT_NEAR(face.GetArea(), 1.0, 1e-12);
		}
		auto CoM = tes.GetCellCM(pt);
		EXPECT_EQ(CoM, mesh[pt]);
		EXPECT_NEAR(tes.GetVolume(pt), 1.0, 1e-12);
//		cout << "Center: " << CoM.x << ", " << CoM.y << ", " << CoM.z << endl;
	}
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