// TestVoroPlusPlus.cpp : Defines the entry point for the console application.
//

#include "gtest/gtest.h"
#include "GeometryCommon/Face.hpp"
#include "GeometryCommon/Vector3D.hpp"
#include "Voronoi/VoroPlusPlus.hpp"

#include <vector>
#include <cstdlib>

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

int main(int argc, char *argv[])
{
	testing::InitGoogleTest(&argc, argv);
	int rc = RUN_ALL_TESTS();

#ifdef _DEBUG
	string s;
	getline(cin, s);
#endif
}