//----------------------------------------------------------------------
// File:			ANN.cpp
// Programmer:		Sunil Arya and David Mount
// Description:		Methods for ANN.h and ANNx.h
// Last modified:	01/27/10 (Version 1.1.2)
//----------------------------------------------------------------------
// Copyright (c) 1997-2010 University of Maryland and Sunil Arya and
// David Mount.  All Rights Reserved.
// 
// This software and related documentation is part of the Approximate
// Nearest Neighbor Library (ANN).  This software is provided under
// the provisions of the Lesser GNU Public License (LGPL).  See the
// file ../ReadMe.txt for further information.
// 
// The University of Maryland (U.M.) and the authors make no
// representations about the suitability or fitness of this software for
// any purpose.  It is provided "as is" without express or implied
// warranty.
//----------------------------------------------------------------------
// History:
//	Revision 0.1  03/04/98
//		Initial release
//	Revision 1.0  04/01/05
//		Added performance counting to annDist()
//	Revision 1.1.2  01/27/10
//		Fixed minor compilation bugs for new versions of gcc
//----------------------------------------------------------------------

#include <cstdlib>						// C standard lib defs
#include "ANNx.h"					// all ANN includes
#include "ANNperf.h"				// ANN performance 

using namespace std;					// make std:: accessible

//----------------------------------------------------------------------
//	Point methods
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	Distance utility.
//		(Note: In the nearest neighbor search, most distances are
//		computed using partial distance calculations, not this
//		procedure.)
//----------------------------------------------------------------------

ANNdist annDist(						// interpoint squared distance
	int					dim,
	ANNpoint const&			p,
	ANNpoint const&			q)
{
	register int d;
	register ANNcoord dist;
	register ANNcoord diff;

	dist = 0;
//#ifdef __INTEL_COMPILER
//#pragma ivdep
//#endif
	for (d = 0; d < dim; d++) 
	{
		diff = p[d] - q[d];
		dist = ANN_SUM(dist, ANN_POW(diff));
	}
	ANN_FLOP(3*dim)					// performance counts
	ANN_PTS(1)
	ANN_COORD(dim)
	return dist;
}

ANNdist annDist2(						// interpoint squared distance
	int					dim,
	ANNpoint	const&		p,
	ANNpoint const&			q,ANNpoint const& cm)
{
	ANNcoord dist=0;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (int d = 0; d < dim; d++)
	{
		ANNcoord diff = std::max(std::abs(p[d] - cm[d]),std::abs(q[d]-cm[d]));
		dist += diff*diff;
	}
	return dist;
}


//----------------------------------------------------------------------
//	annPrintPoint() prints a point to a given output stream.
//----------------------------------------------------------------------

void annPrintPt(						// print a point
	ANNpoint			pt,				// the point
	std::ostream		&out)			// output stream
{
	for (int j = 0; j < 3; j++) {
		out << pt[j];
		if (j < 3-1) out << " ";
	}
}

//----------------------------------------------------------------------
//	Point allocation/deallocation:
//
//		Because points (somewhat like strings in C) are stored
//		as pointers.  Consequently, creating and destroying
//		copies of points may require storage allocation.  These
//		procedures do this.
//
//		annAllocPt() and annDeallocPt() allocate a deallocate
//		storage for a single point, and return a pointer to it.
//
//		annAllocPts() allocates an array of points as well a place
//		to store their coordinates, and initializes the points to
//		point to their respective coordinates.  It allocates point
//		storage in a contiguous block large enough to store all the
//		points.  It performs no initialization.
//
//		annDeallocPts() should only be used on point arrays allocated
//		by annAllocPts since it assumes that points are allocated in
//		a block.
//
//		annCopyPt() copies a point taking care to allocate storage
//		for the new point.
//
//		annAssignRect() assigns the coordinates of one rectangle to
//		another.  The two rectangles must have the same dimension
//		(and it is not possible to test this here).
//----------------------------------------------------------------------
   
ANNpointArray annAllocPts(int n, int /*dim*/)		// allocate n pts in dim
{
	ANNpointArray pa(static_cast<size_t>(n));
	return pa;
}

												// assign one rect to another
void annAssignRect(ANNorthRect &dest, const ANNorthRect &source)
{
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (int i = 0; i < 3; i++) {
		dest.lo[i] = source.lo[i];
		dest.hi[i] = source.hi[i];
	}
}

												// is point inside rectangle?
ANNbool ANNorthRect::inside(ANNpoint p)
{
	for (int i = 0; i < 3; i++) {
		if (p[i] < lo[i] || p[i] > hi[i]) return ANNfalse;
	}
	return ANNtrue;
}

//----------------------------------------------------------------------
//	Error handler
//----------------------------------------------------------------------

void annError(const char* msg, ANNerr level)
{
	if (level == ANNabort) {
		cerr << "ANN: ERROR------->" << msg << "<-------------ERROR\n";
		exit(1);
	}
	else {
		cerr << "ANN: WARNING----->" << msg << "<-------------WARNING\n";
	}
}

//----------------------------------------------------------------------
//	Limit on number of points visited
//		We have an option for terminating the search early if the
//		number of points visited exceeds some threshold.  If the
//		threshold is 0 (its default)  this means there is no limit
//		and the algorithm applies its normal termination condition.
//		This is for applications where there are real time constraints
//		on the running time of the algorithm.
//----------------------------------------------------------------------

int	ANNmaxPtsVisited = 0;	// maximum number of pts visited
int	ANNptsVisited;			// number of pts visited in search

//----------------------------------------------------------------------
//	Global function declarations
//----------------------------------------------------------------------

void annMaxPtsVisit(			// set limit on max. pts to visit in search
	int					maxPts)			// the limit
{
	ANNmaxPtsVisited = maxPts;
}
