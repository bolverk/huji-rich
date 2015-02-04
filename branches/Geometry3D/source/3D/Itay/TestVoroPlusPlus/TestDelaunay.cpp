#include "gtest/gtest.h"
#include "Voronoi/TetGenDelaunay.hpp"
#include "Voronoi/Subcube.hpp"
#include "Utilities/assert.hpp"

#include <unordered_set>
#include <set>

using namespace std;

void BadSubcube()
{
	Subcube s2("--Z");
}

TEST(Subcube, Construction)
{
	Subcube s1("--+");
	ASSERT_THROW(BadSubcube(), invalid_argument);
}

TEST(Subcube, Sets)
{
	set<Subcube> set;
	set.insert(Subcube("---"));
	set.insert(Subcube("--+"));
	set.insert(Subcube("   "));
	set.insert(Subcube("--+"));
	set.insert(Subcube("   "));
	ASSERT_EQ(set.size(), 3);

	unordered_set<Subcube> hash;
	hash.insert(set.begin(), set.end());
	hash.insert(Subcube("+++"));
	hash.insert(Subcube("   "));
	ASSERT_EQ(hash.size(), 4); 
}

OuterBoundary3D InitBoundary()
{
	return OuterBoundary3D(OuterBoundary3D::RECTANGULAR, Vector3D(1, 1, 1), Vector3D(0, 0, 0));
}

TEST(SubCube, DistanceToFace)
{
	OuterBoundary3D boundary = InitBoundary();
	Vector3D pt(0.2, 0.5, 0.9);

	ASSERT_NEAR(distance(pt, boundary, "-  "), 0.2, 1e-8);
	ASSERT_NEAR(distance(pt, boundary, "+  "), 0.8, 1e-8);
	ASSERT_NEAR(distance(pt, boundary, " - "), 0.5, 1e-8);
	ASSERT_NEAR(distance(pt, boundary, " + "), 0.5, 1e-8);
	ASSERT_NEAR(distance(pt, boundary, "  -"), 0.9, 1e-8);
	ASSERT_NEAR(distance(pt, boundary, "  +"), 0.1, 1e-8);
}

TEST(SubCube, DistanceToPoint)
{
	OuterBoundary3D boundary = InitBoundary();
	Vector3D pt(0.2, 0.5, 0.9);

	ASSERT_NEAR(distance(pt, boundary, "---"), 1.0488, 1e-4);
	ASSERT_NEAR(distance(pt, boundary, "+++"), 0.94868, 1e-5);
	ASSERT_NEAR(distance(pt, boundary, "-+-"), 1.0488, 1e-4);  // Same distance because y is 0.5
	ASSERT_NEAR(distance(pt, boundary, "+-+"), 0.94868, 1e-4); // Ditto
	ASSERT_NEAR(distance(pt, boundary, "--+"), 0.54772, 1e-5);
}

TEST(SubCube, DistanceEdge)
{
	OuterBoundary3D boundary = InitBoundary();
	Vector3D pt(0.2, 0.5, 0.9);

	ASSERT_NEAR(distance(pt, boundary, "-- "), 0.53851, 1e-5);
	ASSERT_NEAR(distance(pt, boundary, "+ +"), 0.80623, 1e-5);
	ASSERT_NEAR(distance(pt, boundary, " -+"), 0.50990, 1e-5);
}