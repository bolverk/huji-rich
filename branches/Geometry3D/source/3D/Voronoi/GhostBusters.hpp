//\file GhostBusters.hpp
//\brief Contains the classes that generate and deal with ghost points
//\author Itay Zandbank

#ifndef GHOST_BUSTER_HPP
#define GHOST_BUSTER_HPP 1

#include <set>
#include <vector>
#include <unordered_set>
#include "../GeometryCommon/Vector3D.hpp"
#include "../GeometryCommon/OuterBoundary3D.hpp"
#include "../GeometryCommon/VectorRepository.hpp"
#include "Delaunay.hpp"

//\brief The abstract GhostBuster class. A GhostBuster class is really just a Functor object.
class GhostBuster
{
public:
	virtual std::set<VectorRef> operator()(const Delaunay &del, const OuterBoundary3D &boundary) const = 0;
};

//\brief A simple implementation that copies each point 26 times
class BruteForceGhostBuster : public GhostBuster
{
private:
	static std::set<Subcube> _subcubes;
	static std::set<Subcube> RelevantSubcubes();
public:
	virtual std::set<VectorRef> operator()(const Delaunay &del, const OuterBoundary3D &boundary) const;
};

//\brief An implementation that checks only points that are on edge-faces.
class CloseToBoundaryGhostBuster : public GhostBuster
{
public:
	virtual std::set<VectorRef> operator()(const Delaunay &del, const OuterBoundary3D &boundary) const;

private:
	unordered_set<int> FindOuterTetrahedra(const Delaunay &del) const ;
	unordered_set<int> FindEdgeTetrahedra(const Delaunay &del, const unordered_set<int>& outerTetrahedra) const;

	typedef unordered_map<VectorRef, unordered_set<Subcube>> breach_map;
	breach_map FindHullBreaches(const Delaunay &del, 
		const unordered_set<int>& edgeTetrahedra, 
		const unordered_set<int> &outerTetrahedra, 
		const OuterBoundary3D &boundary) const;

	template<typename T>
	static bool contains(const unordered_set<T> &set, const T &val)
	{
		return set.find(val) != set.end();
	}
};
#endif