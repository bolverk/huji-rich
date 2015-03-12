/*
 \file TetGenDelaunay.hpp
 \brief a TetGen based implementation of the Delaunay abstract class
 \author Itay Zandbank
 */

#ifndef TETGENDELAUNAY_HPP
#define TETGENDELAUNAY_HPP

#include "Delaunay.hpp"
#include "../GeometryCommon/Face.hpp"

class TetGenImpl;
class TetGenDelaunay : public Delaunay
{
public:
	TetGenDelaunay(const std::vector<VectorRef> &points, const Tetrahedron &bigTetrahedron, bool runVoronoi=true);
	
	std::vector<size_t> GetVoronoiCellFaces(size_t cellNum) const;
	std::vector<Face> GetVoronoiFaces() const;

protected:
	void RunDelaunay();
	friend class TetGenImpl;   // This is the PIMPL pattern, although the PIMPL is generated once per call to Run,
							   // so we don't need to hold a pointer to it.

	virtual void FillEdges(); 
	virtual void FillNeighbors();

	std::vector<size_t> OrderNeighbors(const std::vector<size_t> &tetrahedra);

	// Voronoi output
	std::vector<std::vector<int>> _voronoiCellFaces;  // Use ints for all the vectors, because that's what tetgen uses internally.
	std::vector<std::vector<int>> _voronoiFaceEdges;  // Sometimes -1 is used as a marker
	std::vector<std::pair<int, int>> _voronoiFaceNeighbors; // Neighbors of each face
	std::vector<std::pair<int, int>> _voronoiEdges;
	std::vector<Vector3D> _voronoiVertices;

	Face GetVoronoiFace(size_t faceNum) const;
	bool _runVoronoi;
};

#endif // TETGEN_DELAUNAY_HPP