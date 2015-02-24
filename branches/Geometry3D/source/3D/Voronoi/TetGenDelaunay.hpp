/*
 \file TetGenDelaunay.hpp
 \brief a TetGen based implementation of the Delaunay abstract class
 \author Itay Zandbank
 */

#ifndef TETGENDELAUNAY_HPP
#define TETGENDELAUNAY_HPP

#include "Delaunay.hpp"

class TestGenImpl;
class TetGenDelaunay : public Delaunay
{
public:
	TetGenDelaunay(const std::vector<VectorRef> &points, const Tetrahedron &bigTetrahedron);
	
protected:
	void RunDelaunay();
	friend class TetGenImpl;   // This is the PIMPL pattern, although the PIMPL is generated once per call to Run,
							   // so we don't need to hold a pointer to it.

	virtual void FillEdges(); 
	std::vector<int> OrderNeighbors(const std::vector<int> &tetrahedra);

	std::vector<std::vector<int>> _neighbors;  // _neighbors[i] -> the 4 neighbors of tetrahedron i
};

#endif // TETGEN_DELAUNAY_HPP