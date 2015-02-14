#include "gtest/gtest.h"
#include "GeometryCommon/Face.hpp"
#include "GeometryCommon/Vector3D.hpp"
#include "Utilities/assert.hpp"
#include "GeometryCommon/OuterBoundary3D.hpp"
#include "../misc/universal_error.hpp"

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

	Vector3D vec3 = vec2;
	ASSERT_EQ(vec3, vec2);
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
	Face face(vertices);

	ASSERT_EQ(vertices, face.vertices);

	Face face2(vertices);
	face2.AddNeighbor(5, false);
	face2.AddNeighbor(3, true);
	ASSERT_EQ(vertices, face2.vertices);
	ASSERT_EQ(face2.Neighbor1()->GetCell(), 5);
	ASSERT_EQ(face2.Neighbor1()->IsOverlapping(), false);
	ASSERT_EQ(face2.Neighbor2()->GetCell(), 3);
	ASSERT_EQ(face2.Neighbor2()->IsOverlapping(), true);
}

TEST(Geometry3D, Face_Area)
{
	Vector3D v1(0, 0, 0), v2(1, 0, 0), v3(1, 1, 0), v4(0, 1, 0);
	vector<Vector3D> vertices{ v1, v2, v3, v4 };
	Face square(vertices);
	ASSERT_EQ(square.GetArea(), 1);

	Vector3D v5(0, 0, 0), v6(0, 2, 0), v7(2, 0, 0);
	vector<Vector3D> vertices2{ v5, v6, v7 };
	Face triangle(vertices2);
	ASSERT_EQ(triangle.GetArea(), 2);
}

TEST(Geometry3D, Face_Neighbors)
{
	Vector3D v1(0, 0, 0), v2(1, 0, 0), v3(1, 1, 0), v4(0, 1, 0);
	vector<Vector3D> vertices{ v1, v2, v3, v4 };
	Face face1(vertices);

	face1.AddNeighbor(1, true);
	face1.AddNeighbor(2, true);
	face1.AddNeighbor(2, true);
	face1.AddNeighbor(1, true);
	ASSERT_THROW(face1.AddNeighbor(3), UniversalError);
	ASSERT_THROW(face1.AddNeighbor(1, false), UniversalError);
}

TEST(Geometry3D, Face_Identical)
{
	Vector3D v1(0, 0, 0), v2(1, 0, 0), v3(1, 1, 0), v4(0, 1, 0);
	vector<Vector3D> vertices1{ v1, v2, v3, v4 }, vertices2{ v2, v3, v4, v1 }, vertices3{ v4, v3, v2, v1 }, vertices4{ v1, v2, v3 }, vertices5{ v2, v1, v3, v4 };
	Face face1(vertices1), face2(vertices2), face3(vertices3), face4(vertices4);

	ASSERT_TRUE(face1.IdenticalTo(vertices2));
	ASSERT_TRUE(face1.IdenticalTo(vertices3));
	ASSERT_FALSE(face1.IdenticalTo(vertices4));
	ASSERT_TRUE(face2.IdenticalTo(vertices3));
	ASSERT_FALSE(face2.IdenticalTo(vertices4));
	ASSERT_FALSE(face3.IdenticalTo(vertices4));
	ASSERT_FALSE(face2.IdenticalTo(vertices4));
}

TEST(Geometry3D, Face_Order)
{
	Vector3D v1(0, 1, 0), v2(1, 1, 0), v3(1, 0.5, 0), v4(1.2, 0, 0), v5(0.8, -0.8, 0), v6(0, -1, 0);

	vector<Vector3D> orig{ v1, v2, v3, v4, v5, v6};
	vector<Vector3D> reshuffled(orig);
	random_shuffle(reshuffled.begin(), reshuffled.end());

	Face face(reshuffled);
	face.ReorderVertices();

	ASSERT_TRUE(face.IdenticalTo(orig));
}

void CreateFaultyBoundary1()
{
	RectangularBoundary3D(Vector3D(1, 1, 1), Vector3D(0, 2, 1));
}

TEST(Geometry3D, OuterBoundry3D)
{
	RectangularBoundary3D b1(Vector3D(0, 0, 0), Vector3D(-1, -1, -1));
	ASSERT_THROW(CreateFaultyBoundary1(), invalid_argument);
}
