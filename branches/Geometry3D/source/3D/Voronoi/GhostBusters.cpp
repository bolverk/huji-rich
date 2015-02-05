//\file GhostBusters.cpp
//\brief Implementation of the common Ghost Busters
//\author Itay Zandbank

#include "GhostBusters.hpp"
#include <set>
#include <unordered_set>

using namespace std;

set<Vector3D> BruteForceGhostBuster::operator()(const Delaunay &del, const OuterBoundary3D &boundary) const
{
	set<Vector3D> ghosts;
	for (auto vec : del.InputPoints())
		for (auto subcube : Subcube::all())
			ghosts.insert(boundary.ghost(vec, subcube));

	return ghosts;
}

set<Vector3D> CloseToBoundaryGhostBuster::operator()(const Delaunay &del, const OuterBoundary3D &boundary) const
{
	unordered_set<int> outer = FindOuterTetrahedra(del);
	unordered_set<int> edge = FindEdgeTetrahedra(del, outer);
	map<Vector3D, unordered_set<Subcube>> breaches = FindHullBreaches(del, edge, outer, boundary);

	set<Vector3D> ghosts;
	if (breaches.empty())
		return ghosts;

	for (map<Vector3D, unordered_set<Subcube>>::iterator it = breaches.begin(); it != breaches.end(); it++)
	{
		Vector3D pt = it->first;
		for (unordered_set<Subcube>::iterator itSubcube = it->second.begin(); itSubcube != it->second.end(); itSubcube++)
		{
			Vector3D ghost = boundary.ghost(pt, *itSubcube);
			ghosts.insert(ghost);
		}
	}

	return ghosts;
}

CloseToBoundaryGhostBuster::breach_map CloseToBoundaryGhostBuster::FindHullBreaches(const Delaunay &del, const unordered_set<int>& edgeTetrahedra,
	const unordered_set<int> &outerTetrahedra, const OuterBoundary3D &boundary) const
{
	breach_map result;
	for (unordered_set<int>::iterator itTetra = edgeTetrahedra.begin(); itTetra != edgeTetrahedra.end(); itTetra++)
	{
		Tetrahedron t = del[*itTetra];
		for (int iv = 0; iv < 4; iv++)
		{
			Vector3D pt = t[iv];
			if (result.find(pt) != result.end())
				continue;

			vector<int> tetrahedraIndices = del.VertexNeighbors(pt);
			unordered_set<Subcube> breaches;

			for (vector<int>::iterator itTetrahedron = tetrahedraIndices.begin(); itTetrahedron != tetrahedraIndices.end(); itTetrahedron++)
			{
				if (contains(outerTetrahedra, *itTetrahedron))
					continue;

				Tetrahedron tetrahedron = del[*itTetrahedron];
				for (set<Subcube>::iterator itSubcube = Subcube::all().begin(); itSubcube != Subcube::all().end(); itSubcube++)
				{
					if (boundary.distance(tetrahedron.center(), *itSubcube) < tetrahedron.radius())
						breaches.insert(*itSubcube);
				}
			}
			result[pt] = breaches;
		}
	}

	return result;
}

//\brief Find all the Outer Tetrahedra
//\remark An Outer Tetrahedron is a tetrahedron that has a vertex in the big tetrahedron
unordered_set<int> CloseToBoundaryGhostBuster::FindOuterTetrahedra(const Delaunay &del) const
{
	unordered_set<int> result;

	for (int i = 0; i < 4; i++)
	{
		Vector3D pt = del.BigTetrahedron()[i];
		vector<int> tetrahedra = del.VertexNeighbors(pt);
		result.insert(tetrahedra.begin(), tetrahedra.end());
	}

	return result;
}

//\brief Find all the Edge Tetrahedra
//\remark An Edge Tetrahedron has a vertex that belongs to an Outer Tetrahedron. Edge Tetrahedra are not Outer Tetrahedra
unordered_set<int> CloseToBoundaryGhostBuster::FindEdgeTetrahedra(const Delaunay &del, const unordered_set<int>& outerTetrahedra) const
{
	unordered_set<int> edgeTetrahedra;

	// Go over all the outer tetrahedra
	for (unordered_set<int>::const_iterator it = outerTetrahedra.begin(); it != outerTetrahedra.end(); it++)
	{
		Tetrahedron t = del[*it];
		// And all their vertices
		for (int i = 0; i < 4; i++)
		{
			Vector3D pt = t[i];
			if (del.IsBigTetrahedron(pt))  // Ignore the BigTetrahedron - all tetrahedra touching this vertex are Outer and not Edge
				continue;
			// Find all the neighbor tetrahedra of pt
			vector<int> tetrahedra = del.VertexNeighbors(pt);
			// Add them to the set
			for (vector<int>::iterator ptIt = tetrahedra.begin(); ptIt != tetrahedra.end(); ptIt++)
				if (!contains(outerTetrahedra, *ptIt))  // But only if their not Outer
					edgeTetrahedra.insert(*ptIt);
		}
	}

	return edgeTetrahedra;
}
