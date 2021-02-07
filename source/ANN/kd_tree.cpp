//----------------------------------------------------------------------
// File:			kd_tree.cpp
// Programmer:		Sunil Arya and David Mount
// Description:		Basic methods for kd-trees.
// Last modified:	01/04/05 (Version 1.0)
//----------------------------------------------------------------------
// Copyright (c) 1997-2005 University of Maryland and Sunil Arya and
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
//		Increased aspect ratio bound (ANN_AR_TOOBIG) from 100 to 1000.
//		Fixed leaf counts to count trivial leaves.
//		Added optional pa, pi arguments to Skeleton kd_tree constructor
//			for use in load constructor.
//		Added annClose() to eliminate KD_TRIVIAL memory leak.
//----------------------------------------------------------------------

#include <xmmintrin.h>
#include "kd_tree.h"					// kd-tree declarations
#include "kd_split.h"					// kd-tree splitting rules
#include "kd_util.h"					// kd-tree utilities
#include "ANNperf.h"				// performance evaluation
#include <cmath>
#include <cfloat>
#include <algorithm>
//----------------------------------------------------------------------
//	Global data
//
//	For some splitting rules, especially with small bucket sizes,
//	it is possible to generate a large number of empty leaf nodes.
//	To save storage we allocate a single trivial leaf node which
//	contains no points.  For messy coding reasons it is convenient
//	to have it reference a trivial point index.
//
//	KD_TRIVIAL is allocated when the first kd-tree is created.  It
//	must *never* deallocated (since it may be shared by more than
//	one tree).
//----------------------------------------------------------------------
static int				IDX_TRIVIAL[] = { 0 };	// trivial point index
ANNkd_leaf* KD_TRIVIAL = NULL;		// trivial leaf node

//----------------------------------------------------------------------
//	Printing the kd-tree 
//		These routines print a kd-tree in reverse inorder (high then
//		root then low).  (This is so that if you look at the output
//		from the right side it appear from left to right in standard
//		inorder.)  When outputting leaves we output only the point
//		indices rather than the point coordinates. There is an option
//		to print the point coordinates separately.
//
//		The tree printing routine calls the printing routines on the
//		individual nodes of the tree, passing in the level or depth
//		in the tree.  The level in the tree is used to print indentation
//		for readability.
//----------------------------------------------------------------------

namespace
{
	double fastsqrt(double x)
	{
		double res = static_cast<double>(_mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(static_cast<float>(x)))));
		return x * res * (1.5 - 0.5 * res * res * x);
	}

	void CrossProduct(ANNpoint v1, ANNpoint v2, ANNpoint& res)
	{
		res[0] = v1[1] * v2[2] - v1[2] * v2[1];
		res[1] = v1[2] * v2[1] - v1[1] * v2[2];
		res[2] = v1[1] * v2[1] - v1[1] * v2[1];
	}

	double ScalarProduct(ANNpoint v1, ANNpoint v2)
	{
		return v1[0] * v2[0] + v1[1] * v2[1] + v2[2] * v1[2];
	}

	bool PointInBox(ANNorthRect const& bb, ANNpoint const& p)
	{
		int res = 0;
		for (size_t i = 0; i < 3; ++i)
		{
			if (p[i]<bb.lo[i] || p[i] > bb.hi[i])
				++res;
		}
		return res == 0;
	}

	bool PointInface(ANNpointArray face, size_t Nface, ANNpoint point)
	{
		ANNpoint normal, temp0, temp1, res;
		temp0[0] = face[0][0] - point[0];
		temp0[1] = face[0][1] - point[1];
		temp0[2] = face[0][2] - point[2];
		temp1[0] = face[1][0] - point[0];
		temp1[1] = face[1][1] - point[1];
		temp1[2] = face[1][2] - point[2];
		CrossProduct(temp0, temp1, normal);
		--Nface;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
		for (size_t i = 0; i < Nface; i++)
		{
			temp0[0] = face[i + 1][0] - point[0];
			temp0[1] = face[i + 1][1] - point[1];
			temp0[2] = face[i + 1][2] - point[2];
			temp1[0] = face[(i + 2) % (Nface + 1)][0] - point[0];
			temp1[1] = face[(i + 2) % (Nface + 1)][1] - point[1];
			temp1[2] = face[(i + 2) % (Nface + 1)][2] - point[2];
			CrossProduct(temp0, temp1, res);
			if (ScalarProduct(res, normal) < 0)
				return false;
		}
		return true;
	}
}

void ANNkd_split::print(				// print splitting node
	int level,						// depth of node in tree
	ostream& out)					// output stream
{
	child[ANN_HI]->print(level + 1, out);	// print high child
	out << "    ";
	for (int i = 0; i < level; i++)		// print indentation
		out << "..";
	out << "Split cd=" << cut_dim << " cv=" << cut_val;
	out << " lbnd=" << cd_bnds[ANN_LO];
	out << " hbnd=" << cd_bnds[ANN_HI];
	out << "\n";
	child[ANN_LO]->print(level + 1, out);	// print low child
}

void ANNkd_leaf::print(					// print leaf node
	int level,						// depth of node in tree
	ostream& out)					// output stream
{

	out << "    ";
	for (int i = 0; i < level; i++)		// print indentation
		out << "..";

	if (this == KD_TRIVIAL) {			// canonical trivial leaf node
		out << "Leaf (trivial)\n";
	}
	else {
		out << "Leaf n=" << n_pts << " <";
		for (int j = 0; j < n_pts; j++) {
			out << bkt[j];
			if (j < n_pts - 1) out << ",";
		}
		out << ">\n";
	}
}

void ANNkd_tree::Print(					// print entire tree
	ANNbool with_pts,				// print points as well?
	ostream& out)					// output stream
{
	out << "ANN Version " << ANNversion << "\n";
	if (with_pts) {						// print point coordinates
		out << "    Points:\n";
		for (int i = 0; i < n_pts; i++) {
			out << "\t" << i << ": ";
			annPrintPt(pts->operator[](i), out);
			out << "\n";
		}
	}
	if (root == NULL)					// empty tree?
		out << "    Null tree.\n";
	else {
		root->print(0, out);			// invoke printing at root
	}
}

//----------------------------------------------------------------------
//	kd_tree statistics (for performance evaluation)
//		This routine compute various statistics information for
//		a kd-tree.  It is used by the implementors for performance
//		evaluation of the data structure.
//----------------------------------------------------------------------

#define MAX(a,b)		((a) > (b) ? (a) : (b))

void ANNkdStats::merge(const ANNkdStats& st)	// merge stats from child 
{
	n_lf += st.n_lf;			n_tl += st.n_tl;
	n_spl += st.n_spl;			n_shr += st.n_shr;
	depth = MAX(depth, st.depth);
	sum_ar += st.sum_ar;
}

//----------------------------------------------------------------------
//	Update statistics for nodes
//----------------------------------------------------------------------

const double ANN_AR_TOOBIG = 1000;				// too big an aspect ratio

void ANNkd_leaf::getStats(						// get subtree statistics
	ANNkdStats& st,					// stats (modified)
	ANNorthRect& bnd_box)				// bounding box
{
	st.reset();
	st.n_lf = 1;								// count this leaf
	if (this == KD_TRIVIAL) st.n_tl = 1;		// count trivial leaf
	double ar = annAspectRatio(3, bnd_box);	// aspect ratio of leaf
												// incr sum (ignore outliers)
	st.sum_ar += float(ar < ANN_AR_TOOBIG ? ar : ANN_AR_TOOBIG);
}

void ANNkd_split::getStats(						// get subtree statistics
	ANNkdStats& st,					// stats (modified)
	ANNorthRect& bnd_box)				// bounding box
{
	ANNkdStats ch_stats;						// stats for children
												// get stats for low child
	ANNcoord hv = bnd_box.hi[cut_dim];			// save box bounds
	bnd_box.hi[cut_dim] = cut_val;				// upper bound for low child
	ch_stats.reset();							// reset
	child[ANN_LO]->getStats(ch_stats, bnd_box);
	st.merge(ch_stats);							// merge them
	bnd_box.hi[cut_dim] = hv;					// restore bound
												// get stats for high child
	ANNcoord lv = bnd_box.lo[cut_dim];			// save box bounds
	bnd_box.lo[cut_dim] = cut_val;				// lower bound for high child
	ch_stats.reset();							// reset
	child[ANN_HI]->getStats(ch_stats, bnd_box);
	st.merge(ch_stats);							// merge them
	bnd_box.lo[cut_dim] = lv;					// restore bound

	st.depth++;									// increment depth
	st.n_spl++;									// increment number of splits
}

//----------------------------------------------------------------------
//	getStats
//		Collects a number of statistics related to kd_tree or
//		bd_tree.
//----------------------------------------------------------------------

void ANNkd_tree::getStats(						// get tree statistics
	ANNkdStats& st)					// stats (modified)
{
	st.reset(3, n_pts, bkt_size);				// reset stats
												// create bounding box
	ANNorthRect bnd_box(bnd_box_lo, bnd_box_hi);
	if (root != NULL) {							// if nonempty tree
		root->getStats(st, bnd_box);		// get statistics
		st.avg_ar = st.sum_ar / (float)st.n_lf;		// average leaf asp ratio
	}
}

//----------------------------------------------------------------------
//	kd_tree destructor
//		The destructor just frees the various elements that were
//		allocated in the construction process.
//----------------------------------------------------------------------

ANNkd_tree::~ANNkd_tree()				// tree destructor
{
	if (root != NULL) delete root;
	if (pidx != NULL) delete[] pidx;
}

//----------------------------------------------------------------------
//	This is called with all use of ANN is finished.  It eliminates the
//	minor memory leak caused by the allocation of KD_TRIVIAL.
//----------------------------------------------------------------------
void annClose()				// close use of ANN
{
	if (KD_TRIVIAL != NULL) {
		delete KD_TRIVIAL;
		KD_TRIVIAL = NULL;
	}
}

//----------------------------------------------------------------------
//	kd_tree constructors
//		There is a skeleton kd-tree constructor which sets up a
//		trivial empty tree.	 The last optional argument allows
//		the routine to be passed a point index array which is
//		assumed to be of the proper size (n).  Otherwise, one is
//		allocated and initialized to the identity.	Warning: In
//		either case the destructor will deallocate this array.
//
//		As a kludge, we need to allocate KD_TRIVIAL if one has not
//		already been allocated.	 (This is because I'm too dumb to
//		figure out how to cause a pointer to be allocated at load
//		time.)
//----------------------------------------------------------------------

void ANNkd_tree::SkeletonTree(			// construct skeleton tree
	int n,							// number of points
	int bs,							// bucket size
	ANNpointArray pa,				// point array
	ANNidxArray pi)					// point indices
{
	n_pts = n;
	bkt_size = bs;
	pts = &pa;							// initialize points array

	root = NULL;						// no associated tree yet

	if (pi == NULL) {					// point indices provided?
		pidx = new ANNidx[n];			// no, allocate space for point indices
		for (int i = 0; i < n; i++) {
			pidx[i] = i;				// initially identity
		}
	}
	else {
		pidx = pi;						// yes, use them
	}

	if (KD_TRIVIAL == NULL)				// no trivial leaf node yet?
		KD_TRIVIAL = new ANNkd_leaf(0, IDX_TRIVIAL);	// allocate it
}

ANNkd_tree::ANNkd_tree(					// basic constructor
	int n,							// number of points
	int bs)							// bucket size
	:n_pts(0), bkt_size(0), pts(0), pidx(0), root(0), bnd_box_lo(ANNpoint()), bnd_box_hi(ANNpoint())
{
	SkeletonTree(n, bs);
}			// construct skeleton tree

//----------------------------------------------------------------------
//	rkd_tree - recursive procedure to build a kd-tree
//
//		Builds a kd-tree for points in pa as indexed through the
//		array pidx[0..n-1] (typically a subarray of the array used in
//		the top-level call).  This routine permutes the array pidx,
//		but does not alter pa[].
//
//		The construction is based on a standard algorithm for constructing
//		the kd-tree (see Friedman, Bentley, and Finkel, ``An algorithm for
//		finding best matches in logarithmic expected time,'' ACM Transactions
//		on Mathematical Software, 3(3):209-226, 1977).  The procedure
//		operates by a simple divide-and-conquer strategy, which determines
//		an appropriate orthogonal cutting plane (see below), and splits
//		the points.  When the number of points falls below the bucket size,
//		we simply store the points in a leaf node's bucket.
//
//		One of the arguments is a pointer to a splitting routine,
//		whose prototype is:
//		
//				void split(
//						ANNpointArray pa,  // complete point array
//						ANNidxArray pidx,  // point array (permuted on return)
//						ANNorthRect &bnds, // bounds of current cell
//						int n,			   // number of points
//						int dim,		   // dimension of space
//						int &cut_dim,	   // cutting dimension
//						ANNcoord &cut_val, // cutting value
//						int &n_lo)		   // no. of points on low side of cut
//
//		This procedure selects a cutting dimension and cutting value,
//		partitions pa about these values, and returns the number of
//		points on the low side of the cut.
//----------------------------------------------------------------------

ANNkd_ptr rkd_tree(				// recursive construction of kd-tree
	ANNpointArray const& pa,				// point array
	ANNidxArray			pidx,			// point indices to store in subtree
	int					n,				// number of points
	int					dim,			// dimension of space
	int					bsp,			// bucket space
	ANNorthRect& bnd_box,		// bounding box for current node
	ANNkd_splitter		splitter)		// splitting routine
{
	if (n <= bsp) {						// n small, make a leaf node
		if (n == 0)						// empty leaf node
			return KD_TRIVIAL;			// return (canonical) empty leaf
		else							// construct the node and return
			return new ANNkd_leaf(n, pidx);
	}
	else {								// n large, make a splitting node
		int cd;							// cutting dimension
		ANNcoord cv;					// cutting value
		int n_lo;						// number on low side of cut
		ANNkd_node* lo, * hi;			// low and high children

										// invoke splitting procedure
		(*splitter)(pa, pidx, bnd_box, n, dim, cd, cv, n_lo);

		ANNcoord lv = bnd_box.lo[cd];	// save bounds for cutting dimension
		ANNcoord hv = bnd_box.hi[cd];

		bnd_box.hi[cd] = cv;			// modify bounds for left subtree
		lo = rkd_tree(					// build left subtree
			pa, pidx, n_lo,			// ...from pidx[0..n_lo-1]
			dim, bsp, bnd_box, splitter);
		bnd_box.hi[cd] = hv;			// restore bounds

		bnd_box.lo[cd] = cv;			// modify bounds for right subtree
		hi = rkd_tree(					// build right subtree
			pa, pidx + n_lo, n - n_lo,// ...from pidx[n_lo..n-1]
			dim, bsp, bnd_box, splitter);
		bnd_box.lo[cd] = lv;			// restore bounds

										// create the splitting node
		ANNkd_split* ptr = new ANNkd_split(cd, cv, lv, hv, lo, hi);

		return ptr;						// return pointer to this node
	}
}

ANNkd_ptr rkd_tree(				// recursive construction of kd-tree
	ANNpointArray const& pa,				// point array
	ANNidxArray			pidx,			// point indices to store in subtree
	vector<double> const& masses,
	std::vector<std::array<double, 6> > const& Qs,
	int					n,				// number of points
	int					dim,			// dimension of space
	int					bsp,			// bucket space
	ANNorthRect& bnd_box,		// bounding box for current node
	ANNkd_splitter		splitter)		// splitting routine
{
	if (n <= bsp) {						// n small, make a leaf node
		if (n == 0)						// empty leaf node
			return KD_TRIVIAL;			// return (canonical) empty leaf
		else							// construct the node and return
			return new ANNkd_leaf(n, pidx, masses[pidx[0]], Qs[pidx[0]], pa[pidx[0]]);
	}
	else {								// n large, make a splitting node
		int cd;							// cutting dimension
		ANNcoord cv;					// cutting value
		int n_lo;						// number on low side of cut
		ANNkd_node* lo, * hi;			// low and high children

										// invoke splitting procedure
		(*splitter)(pa, pidx, bnd_box, n, dim, cd, cv, n_lo);

		ANNcoord lv = bnd_box.lo[cd];	// save bounds for cutting dimension
		ANNcoord hv = bnd_box.hi[cd];

		bnd_box.hi[cd] = cv;			// modify bounds for left subtree
		lo = rkd_tree(					// build left subtree
			pa, pidx, masses, Qs, n_lo,			// ...from pidx[0..n_lo-1]
			dim, bsp, bnd_box, splitter);
		bnd_box.hi[cd] = hv;			// restore bounds

		bnd_box.lo[cd] = cv;			// modify bounds for right subtree
		hi = rkd_tree(					// build right subtree
			pa, pidx + n_lo, masses, Qs, n - n_lo,// ...from pidx[n_lo..n-1]
			dim, bsp, bnd_box, splitter);
		bnd_box.lo[cd] = lv;			// restore bounds

										// create the splitting node
		ANNkd_split* ptr = new ANNkd_split(cd, cv, lv, hv, lo, hi);

		ptr->mass = lo->mass + hi->mass;
		for (size_t i = 0; i < 3; ++i)
			ptr->CM[i] = (lo->CM[i] * lo->mass + hi->CM[i] * hi->mass) / std::max(DBL_MIN * 1e10, ptr->mass);
		double qx = lo->CM[0] - ptr->CM[0];
		double qy = lo->CM[1] - ptr->CM[1];
		double qz = lo->CM[2] - ptr->CM[2];
		double qr2 = qx * qx + qy * qy + qz * qz;
		ptr->Q[0] = lo->Q[0] + lo->mass * (3 * qx * qx - qr2);
		ptr->Q[1] = lo->Q[1] + 3 * lo->mass * qx * qy;
		ptr->Q[2] = lo->Q[2] + 3 * lo->mass * qx * qz;
		ptr->Q[3] = lo->Q[3] + lo->mass * (3 * qy * qy - qr2);
		ptr->Q[4] = lo->Q[4] + 3 * lo->mass * qz * qy;
		qx = hi->CM[0] - ptr->CM[0];
		qy = hi->CM[1] - ptr->CM[1];
		qz = hi->CM[2] - ptr->CM[2];
		qr2 = qx * qx + qy * qy + qz * qz;
		ptr->Q[0] += hi->Q[0] + hi->mass * (3 * qx * qx - qr2);
		ptr->Q[1] += hi->Q[1] + 3 * hi->mass * qx * qy;
		ptr->Q[2] += hi->Q[2] + 3 * hi->mass * qx * qz;
		ptr->Q[3] += hi->Q[3] + hi->mass * (3 * qy * qy - qr2);
		ptr->Q[4] += hi->Q[4] + 3 * hi->mass * qz * qy;
		ptr->Q[5] = -ptr->Q[0] - ptr->Q[3];
		return ptr;						// return pointer to this node
	}
}

ANNkd_ptr rkd_tree(ANNpointArray const& pa, ANNidxArray pidx, vector<double> const& masses, int n, int dim, int bsp, ANNorthRect& bnd_box, ANNkd_splitter splitter)
{
	if (n <= bsp) {						// n small, make a leaf node
		if (n == 0)						// empty leaf node
			return KD_TRIVIAL;			// return (canonical) empty leaf
		else							// construct the node and return
			return new ANNkd_leaf(n, pidx, masses[pidx[0]], pa[pidx[0]]);
	}
	else
	{								// n large, make a splitting node
		int cd;							// cutting dimension
		ANNcoord cv;					// cutting value
		int n_lo;						// number on low side of cut
		ANNkd_node* lo, * hi;			// low and high children

										// invoke splitting procedure
		(*splitter)(pa, pidx, bnd_box, n, dim, cd, cv, n_lo);

		ANNcoord lv = bnd_box.lo[cd];	// save bounds for cutting dimension
		ANNcoord hv = bnd_box.hi[cd];

		bnd_box.hi[cd] = cv;			// modify bounds for left subtree
		lo = rkd_tree(					// build left subtree
			pa, pidx, masses, n_lo,			// ...from pidx[0..n_lo-1]
			dim, bsp, bnd_box, splitter);
		bnd_box.hi[cd] = hv;			// restore bounds

		bnd_box.lo[cd] = cv;			// modify bounds for right subtree
		hi = rkd_tree(					// build right subtree
			pa, pidx + n_lo, masses, n - n_lo,// ...from pidx[n_lo..n-1]
			dim, bsp, bnd_box, splitter);
		bnd_box.lo[cd] = lv;			// restore bounds

										// create the splitting node
		ANNkd_split* ptr = new ANNkd_split(cd, cv, lv, hv, lo, hi);

		ptr->mass = lo->mass + hi->mass;
		for (size_t i = 0; i < 3; ++i)
			ptr->CM[i] = (lo->CM[i] * lo->mass + hi->CM[i] * hi->mass) / std::max(DBL_MIN * 1e10, ptr->mass);
		return ptr;						// return pointer to this node
	}
}


//----------------------------------------------------------------------
// kd-tree constructor
//		This is the main constructor for kd-trees given a set of points.
//		It first builds a skeleton tree, then computes the bounding box
//		of the data points, and then invokes rkd_tree() to actually
//		build the tree, passing it the appropriate splitting routine.
//----------------------------------------------------------------------

ANNkd_tree::ANNkd_tree(					// construct from point array
	ANNpointArray const& pa,				// point array (with at least n pts)
	int					n,				// number of points
	int					bs,				// bucket size
	ANNsplitRule		split)			// splitting method
	:n_pts(0), bkt_size(0), pts(0), pidx(0), root(0), bnd_box_lo(ANNpoint()), bnd_box_hi(ANNpoint())
{
	SkeletonTree(n, bs);			// set up the basic stuff
	pts = &pa;							// where the points are
	if (n == 0) return;					// no points--no sweat

	ANNorthRect bnd_box;			// bounding box for points
	annEnclRect(pa, pidx, n, 3, bnd_box);// construct bounding rectangle
										// copy to tree structure
	bnd_box_lo = bnd_box.lo;
	bnd_box_hi = bnd_box.hi;

	switch (split) {					// build by rule
	case ANN_KD_STD:					// standard kd-splitting rule
		root = rkd_tree(pa, pidx, n, 3, bs, bnd_box, kd_split);
		break;
	case ANN_KD_MIDPT:					// midpoint split
		root = rkd_tree(pa, pidx, n, 3, bs, bnd_box, midpt_split);
		break;
	case ANN_KD_FAIR:					// fair split
		root = rkd_tree(pa, pidx, n, 3, bs, bnd_box, fair_split);
		break;
	case ANN_KD_SUGGEST:				// best (in our opinion)
	case ANN_KD_SL_MIDPT:				// sliding midpoint split
		root = rkd_tree(pa, pidx, n, 3, bs, bnd_box, sl_midpt_split);
		break;
	case ANN_KD_SL_FAIR:				// sliding fair split
		root = rkd_tree(pa, pidx, n, 3, bs, bnd_box, sl_fair_split);
		break;
	default:
		annError("Illegal splitting method", ANNabort);
	}
}

ANNkd_tree::ANNkd_tree(					// construct from point array
	ANNpointArray const& pa,				// point array (with at least n pts)
	vector<double> const& masses,
	std::vector<std::array<double, 6> > const& Qs,
	int					n,				// number of points
	int					bs,				// bucket size
	ANNsplitRule		split)			// splitting method
	:n_pts(0), bkt_size(0), pts(0), pidx(0), root(0), bnd_box_lo(ANNpoint()), bnd_box_hi(ANNpoint())
{
	SkeletonTree(n, bs);			// set up the basic stuff
	pts = &pa;							// where the points are
	if (n == 0) return;					// no points--no sweat

	ANNorthRect bnd_box;			// bounding box for points
	annEnclRect(pa, pidx, n, 3, bnd_box);// construct bounding rectangle
										  // copy to tree structure
	bnd_box_lo = bnd_box.lo;
	bnd_box_hi = bnd_box.hi;

	switch (split) {					// build by rule
	case ANN_KD_STD:					// standard kd-splitting rule
		root = rkd_tree(pa, pidx, masses, Qs, n, 3, bs, bnd_box, kd_split);
		break;
	case ANN_KD_MIDPT:					// midpoint split
		root = rkd_tree(pa, pidx, masses, Qs, n, 3, bs, bnd_box, midpt_split);
		break;
	case ANN_KD_FAIR:					// fair split
		root = rkd_tree(pa, pidx, masses, Qs, n, 3, bs, bnd_box, fair_split);
		break;
	case ANN_KD_SUGGEST:				// best (in our opinion)
	case ANN_KD_SL_MIDPT:				// sliding midpoint split
		root = rkd_tree(pa, pidx, masses, Qs, n, 3, bs, bnd_box, sl_midpt_split);
		break;
	case ANN_KD_SL_FAIR:				// sliding fair split
		root = rkd_tree(pa, pidx, masses, Qs, n, 3, bs, bnd_box, sl_fair_split);
		break;
	default:
		annError("Illegal splitting method", ANNabort);
	}
}

namespace
{
	std::pair<bool, double> ZDistanceToFace(ANNorthRect const& face, ANNorthRect const& bb)
	{
		if (!(face.lo[0] >= bb.hi[0] || face.hi[0] <= bb.lo[0] || face.hi[1] <= bb.lo[1] || face.lo[1] >= bb.hi[1]))
			return std::pair<bool, double>(true, std::min(std::abs(bb.lo[2] - face.hi[2]),
				std::abs(bb.hi[2] - face.lo[2])));
		else
			return std::pair<bool, double>(false, 0);
	}

	std::pair<bool, double> ZDistanceToFaces(std::vector<ANNorthRect> const& faces, ANNorthRect const& bb)
	{
		double min_1 = -1;
		bool intersect = false;
		size_t Nfaces = faces.size();
		for (size_t i = 0; i < Nfaces; ++i)
		{
			std::pair<bool, double> res = ZDistanceToFace(faces[i], bb);
			if (res.first)
			{
				intersect = true;
				min_1 = std::max(min_1, 1.0 / res.second);
			}
		}
		return std::pair<bool, double>(intersect, 1.0 / min_1);
	}

	double DistanceToFace(ANNpointArray const& face, size_t Nface, std::array<double, 3> const& qpoint, double maxdist,
		ANNpoint normal)
	{
		double d = (qpoint[0] - face[0][0]) * normal[0] + (qpoint[1] - face[0][1]) * normal[1] + (qpoint[2] - face[0][2]) * normal[2];
		if (std::abs(d) > maxdist)
			return std::abs(d);
		ANNpoint ptemp;
		ptemp[0] = qpoint[0] - d * normal[0];
		ptemp[1] = qpoint[1] - d * normal[1];
		ptemp[2] = qpoint[2] - d * normal[2];
		if (PointInface(face, Nface, ptemp))
		{
			return std::abs(d);
		}
		double min_face_dist = 0;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
		for (size_t i = 0; i < 3; ++i)
			min_face_dist += (qpoint[i] - face[0][i]) * (qpoint[i] - face[0][i]);
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
		for (size_t j = 1; j < Nface; ++j)
		{
			double temp = 0;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
			for (size_t i = 0; i < 3; ++i)
				temp += (qpoint[i] - face[j][i]) * (qpoint[i] - face[j][i]);
			min_face_dist = std::min(min_face_dist, temp);
		}
		return fastsqrt(min_face_dist);
	}

	double DistanceToFaces(std::vector<ANNpointArray> const& faces, std::vector<size_t>const& Nface, std::array<double, 3> const&
		qpoint, double maxdist, std::vector<ANNpoint>const& normals)
	{
		double res = DistanceToFace(faces[0], Nface[0], qpoint, maxdist, normals[0]);
		size_t Nfaces = faces.size();
		for (size_t i = 1; i < Nfaces; ++i)
		{
			double temp = DistanceToFace(faces[i], Nface[i], qpoint, maxdist, normals[i]);
			res = std::min(res, temp);
		}
		return res;
	}
}



ANNkd_tree::ANNkd_tree(ANNpointArray const& pa, std::vector<double> const& masses, int n, int bs, ANNsplitRule split)
	:n_pts(0), bkt_size(0), pts(0), pidx(0), root(0), bnd_box_lo(ANNpoint()), bnd_box_hi(ANNpoint())
{
	SkeletonTree(n, bs);			// set up the basic stuff
	pts = &pa;							// where the points are
	if (n == 0) return;					// no points--no sweat

	ANNorthRect bnd_box;			// bounding box for points
	annEnclRect(pa, pidx, n, 3, bnd_box);// construct bounding rectangle
										  // copy to tree structure
	bnd_box_lo = bnd_box.lo;
	bnd_box_hi = bnd_box.hi;

	switch (split) {					// build by rule
	case ANN_KD_STD:					// standard kd-splitting rule
		root = rkd_tree(pa, pidx, masses, n, 3, bs, bnd_box, kd_split);
		break;
	case ANN_KD_MIDPT:					// midpoint split
		root = rkd_tree(pa, pidx, masses, n, 3, bs, bnd_box, midpt_split);
		break;
	case ANN_KD_FAIR:					// fair split
		root = rkd_tree(pa, pidx, masses, n, 3, bs, bnd_box, fair_split);
		break;
	case ANN_KD_SUGGEST:				// best (in our opinion)
	case ANN_KD_SL_MIDPT:				// sliding midpoint split
		root = rkd_tree(pa, pidx, masses, n, 3, bs, bnd_box, sl_midpt_split);
		break;
	case ANN_KD_SL_FAIR:				// sliding fair split
		root = rkd_tree(pa, pidx, masses, n, 3, bs, bnd_box, sl_fair_split);
		break;
	default:
		annError("Illegal splitting method", ANNabort);
	}
}

void ANNkd_tree::GetAcc(ANNpoint const& qpoint, ANNpoint& res, double angle2) const
{
	ANNorthRect bb(bnd_box_lo, bnd_box_hi);
	root->GetAcc(qpoint, res, angle2, bb);
}

void  ANNkd_tree::GetAcc(std::vector<ANNpoint>& qpoint,
	std::vector<ANNpoint>& res, double angle2) const
{
	ANNorthRect bb(bnd_box_lo, bnd_box_hi);
	std::array<double, 3> qMin, qMax;
	qMax[0] = qpoint[0][0];
	qMax[1] = qpoint[0][1];
	qMax[2] = qpoint[0][2];
	qMin = qMax;
	size_t N = qpoint.size();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t i = 1; i < N; ++i)
	{
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
		for (size_t j = 0; j < 3; ++j)
		{
			qMax[j] = std::max(qMax[j], qpoint[i][j]);
			qMin[j] = std::min(qMin[j], qpoint[i][j]);
		}
	}
	std::array<double, 4> qCM;
	qCM[0] = 0.5 * (qMax[0] + qMin[0]);
	qCM[1] = 0.5 * (qMax[1] + qMin[1]);
	qCM[2] = 0.5 * (qMax[2] + qMin[2]);
	if (N > 1)
		qCM[3] = fastsqrt((qMax[0] - qMin[0]) * (qMax[0] - qMin[0]) + (qMax[1] - qMin[1]) * (qMax[1] - qMin[1]) +
		(qMax[2] - qMin[2]) * (qMax[2] - qMin[2]));
	else
		qCM[3] = 0;
	root->GetAcc(qpoint, res, angle2, bb, qCM);
	}

void ANNkd_split::GetAcc(std::vector<ANNpoint>& qpoint, std::vector<ANNpoint>& res, double angle2, ANNorthRect& bb,
	std::array<double, 4> const& qCM) const
{
	double lv = bb.lo[cut_dim];
	double hv = bb.hi[cut_dim];
	double maxbox = annDist(3, bb.lo, bb.hi);
	int counter = 0;
	size_t N = qpoint.size();
	double dist_toq = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:dist_toq)
#endif
	for (int i = 0; i < 3; ++i)
		dist_toq += (qCM[i] - CM[i]) * (qCM[i] - CM[i]);
	if (N > 1)
		dist_toq -= 2 * qCM[3] * fastsqrt(dist_toq) - qCM[3] * qCM[3];
	if (dist_toq * angle2 > maxbox)
		counter = 1;
	if (counter > 0)
	{
		for (size_t k = 0; k < N; ++k)
		{
			dist_toq = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:dist_toq)
#endif
			for (int i = 0; i < 3; ++i)
				dist_toq += (qpoint[k][i] - CM[i]) * (qpoint[k][i] - CM[i]);
			if (dist_toq * angle2 > maxbox)
			{
				double r3 = 1.0 / (dist_toq * fastsqrt(dist_toq));
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
				for (int i = 0; i < 3; ++i)
					res[k][i] -= mass * (qpoint[k][i] - CM[i]) * r3;
				double Qfactor = r3 / dist_toq;
				double dx = qpoint[k][0] - CM[0];
				double dy = qpoint[k][1] - CM[1];
				double dz = qpoint[k][2] - CM[2];
				res[k][0] += Qfactor * (dx * Q[0] + dy * Q[1] + dz * Q[2]);
				res[k][1] += Qfactor * (dx * Q[1] + dy * Q[3] + dz * Q[4]);
				res[k][2] += Qfactor * (dx * Q[2] + dy * Q[4] + dz * Q[5]);
				double mrr = dx * dx * Q[0] + dy * dy * Q[3] + dz * dz * Q[5] + 2 * dx * dy * Q[1] + 2 * dx * dz * Q[2] + 2 * dy * dz * Q[4];
				Qfactor *= -5 * mrr / (2 * dist_toq);
				res[k][0] += Qfactor * dx;
				res[k][1] += Qfactor * dy;
				res[k][2] += Qfactor * dz;
			}
			else
			{
				if (child[1] != KD_TRIVIAL)
				{
					bb.lo[cut_dim] = cut_val;
					child[1]->GetAcc(qpoint[k], res[k], angle2, bb);
					bb.lo[cut_dim] = lv;
				}
				if (child[0] != KD_TRIVIAL)
				{
					bb.hi[cut_dim] = cut_val;
					child[0]->GetAcc(qpoint[k], res[k], angle2, bb);
					bb.hi[cut_dim] = hv;
				}
			}
		}
		return;
			}
	if (child[1] != KD_TRIVIAL)
	{
		bb.lo[cut_dim] = cut_val;
		child[1]->GetAcc(qpoint, res, angle2, bb);
		bb.lo[cut_dim] = lv;
	}
	if (child[0] != KD_TRIVIAL)
	{
		bb.hi[cut_dim] = cut_val;
		child[0]->GetAcc(qpoint, res, angle2, bb);
		bb.hi[cut_dim] = hv;
	}
	}

void  ANNkd_split::GetAcc(std::vector<ANNpoint>& qpoint, std::vector<ANNpoint>& res, double angle2, ANNorthRect& bb) const
{
	double lv = bb.lo[cut_dim];
	double hv = bb.hi[cut_dim];
	double maxbox = annDist(3, bb.lo, bb.hi);
	size_t counter = 0;
	size_t N = qpoint.size();
	for (size_t k = 0; k < N; ++k)
	{
		double dist_toq = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:dist_toq)
#endif
		for (int i = 0; i < 3; ++i)
			dist_toq += (qpoint[k][i] - CM[i]) * (qpoint[k][i] - CM[i]);
		if (dist_toq * angle2 > maxbox)
		{
			++counter;
			break;
		}
	}
	if (counter > 0)
	{
		for (size_t k = 0; k < N; ++k)
		{
			double dist_toq = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:dist_toq)
#endif
			for (int i = 0; i < 3; ++i)
				dist_toq += (qpoint[k][i] - CM[i]) * (qpoint[k][i] - CM[i]);
			if (dist_toq * angle2 > maxbox)
			{
				double r3 = 1.0 / (dist_toq * fastsqrt(dist_toq));
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
				for (int i = 0; i < 3; ++i)
					res[k][i] -= mass * (qpoint[k][i] - CM[i]) * r3;
				double Qfactor = r3 / dist_toq;
				double dx = qpoint[k][0] - CM[0];
				double dy = qpoint[k][1] - CM[1];
				double dz = qpoint[k][2] - CM[2];
				res[k][0] += Qfactor * (dx * Q[0] + dy * Q[1] + dz * Q[2]);
				res[k][1] += Qfactor * (dx * Q[1] + dy * Q[3] + dz * Q[4]);
				res[k][2] += Qfactor * (dx * Q[2] + dy * Q[4] + dz * Q[5]);
				double mrr = dx * dx * Q[0] + dy * dy * Q[3] + dz * dz * Q[5] + 2 * dx * dy * Q[1] + 2 * dx * dz * Q[2] + 2 * dy * dz * Q[4];
				Qfactor *= -5 * mrr / (2 * dist_toq);
				res[k][0] += Qfactor * dx;
				res[k][1] += Qfactor * dy;
				res[k][2] += Qfactor * dz;
			}
			else
			{
				if (child[1] != KD_TRIVIAL)
				{
					bb.lo[cut_dim] = cut_val;
					child[1]->GetAcc(qpoint[k], res[k], angle2, bb);
					bb.lo[cut_dim] = lv;
				}
				if (child[0] != KD_TRIVIAL)
				{
					bb.hi[cut_dim] = cut_val;
					child[0]->GetAcc(qpoint[k], res[k], angle2, bb);
					bb.hi[cut_dim] = hv;
				}
			}
		}
		return;
	}
	if (child[1] != KD_TRIVIAL)
	{
		bb.lo[cut_dim] = cut_val;
		child[1]->GetAcc(qpoint, res, angle2, bb);
		bb.lo[cut_dim] = lv;
	}
	if (child[0] != KD_TRIVIAL)
	{
		bb.hi[cut_dim] = cut_val;
		child[0]->GetAcc(qpoint, res, angle2, bb);
		bb.hi[cut_dim] = hv;
	}
}

void ANNkd_split::GetToSendOpticalDepth(std::vector<ANNorthRect> const& faces, vector<ANNkd_ptr>& nodes, double angle2, ANNorthRect& bb)
{
	double maxbox = annDist(3, bb.lo, bb.hi);
	std::pair<bool, double> dist = ZDistanceToFaces(faces, bb);
	if (!dist.first || dist.second * dist.second * angle2 > maxbox)
	{
		nodes.push_back(this);
		return;
	}

	double lv = bb.lo[cut_dim];
	double hv = bb.hi[cut_dim];

	if (child[1] != KD_TRIVIAL)
	{
		bb.lo[cut_dim] = cut_val;
		child[1]->GetToSendOpticalDepth(faces, nodes, angle2, bb);
		bb.lo[cut_dim] = lv;
	}
	if (child[0] != KD_TRIVIAL)
	{
		bb.hi[cut_dim] = cut_val;
		child[0]->GetToSendOpticalDepth(faces, nodes, angle2, bb);
		bb.hi[cut_dim] = hv;
	}
}

void ANNkd_split::GetOpticalDepth(ANNpoint const& qpoint, std::vector<std::pair<double, double> >& res, double angle2, ANNorthRect& bb) const
{
	double lv = bb.lo[cut_dim];
	double hv = bb.hi[cut_dim];
	double maxbox = annDist(3, bb.lo, bb.hi);
	double dist_toq = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:dist_toq)
#endif
	for (int i = 0; i < 3; ++i)
		dist_toq += (qpoint[i] - CM[i]) * (qpoint[i] - CM[i]);
	if (dist_toq * angle2 > maxbox && !PointInBox(bb, qpoint))
	{
		if (qpoint[0] >= bb.lo[0] && qpoint[0] <= bb.hi[0] &&
			qpoint[1] >= bb.lo[1] && qpoint[1] <= bb.hi[1])
		{
			res.push_back(std::pair<double, double>(CM[2], mass / ((bb.hi[0] - bb.lo[0]) *
				(bb.hi[1] - bb.lo[1]))));
		}
		return;
	}
	else
	{
		if (child[1] != KD_TRIVIAL)
		{
			bb.lo[cut_dim] = cut_val;
			child[1]->GetOpticalDepth(qpoint, res, angle2, bb);
			bb.lo[cut_dim] = lv;
		}
		if (child[0] != KD_TRIVIAL)
		{
			bb.hi[cut_dim] = cut_val;
			child[0]->GetOpticalDepth(qpoint, res, angle2, bb);
			bb.hi[cut_dim] = hv;
		}
	}
}

void ANNkd_split::GetAcc(ANNpoint const& qpoint, ANNpoint& res, double angle2, ANNorthRect& bb) const
{
	double maxbox = annDist(3, bb.lo, bb.hi);
	double dist_toq = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:dist_toq)
#endif
	for (int i = 0; i < 3; ++i)
		dist_toq += (qpoint[i] - CM[i]) * (qpoint[i] - CM[i]);
	if (dist_toq * angle2 > maxbox)
	{
		double r3 = 1.0 / (dist_toq * fastsqrt(dist_toq));
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
		for (int i = 0; i < 3; ++i)
			res[i] -= mass * (qpoint[i] - CM[i]) * r3;
		const double Qfactor = r3 / dist_toq;
		const double dx = qpoint[0] - CM[0];
		const double dy = qpoint[1] - CM[1];
		const double dz = qpoint[2] - CM[2];
		res[0] += Qfactor * (dx * Q[0] + dy * Q[1] + dz * Q[2]);
		res[1] += Qfactor * (dx * Q[1] + dy * Q[3] + dz * Q[4]);
		res[2] += Qfactor * (dx * Q[2] + dy * Q[4] + dz * Q[5]);
		double mrr = dx * dx * Q[0] + dy * dy * Q[3] + dz * dz * Q[5] + 2 * dx * dy * Q[1] + 2 * dx * dz * Q[2] + 2 * dy * dz * Q[4];
		const double Qfactor2 = -5 * mrr * Qfactor / (2 * dist_toq);
		res[0] += Qfactor2 * dx;
		res[1] += Qfactor2 * dy;
		res[2] += Qfactor2 * dz;
		return;
	}

	double lv = bb.lo[cut_dim];
	double hv = bb.hi[cut_dim];

	if (child[1] != KD_TRIVIAL)
	{
		bb.lo[cut_dim] = cut_val;
		child[1]->GetAcc(qpoint, res, angle2, bb);
		bb.lo[cut_dim] = lv;
	}
	if (child[0] != KD_TRIVIAL)
	{
		bb.hi[cut_dim] = cut_val;
		child[0]->GetAcc(qpoint, res, angle2, bb);
		bb.hi[cut_dim] = hv;
	}
}

void ANNkd_leaf::GetToSendOpticalDepth(std::vector<ANNorthRect> const& /*faces*/, vector<ANNkd_ptr>& nodes, double /*angle2*/, ANNorthRect& /*bb*/)
{
	nodes.push_back(this);
}

void ANNkd_leaf::GetOpticalDepth(ANNpoint const& qpoint, std::vector<std::pair<double, double> >& res, double /*angle2*/, ANNorthRect& bb) const
{
	if (qpoint[0] >= bb.lo[0] && qpoint[0] <= bb.hi[0] &&
		qpoint[1] >= bb.lo[1] && qpoint[1] <= bb.hi[1])
	{
		res.push_back(std::pair<double, double>(CM[2], mass / ((bb.hi[0] - bb.lo[0]) *
			(bb.hi[1] - bb.lo[1]))));
	}
	return;
}

void ANNkd_leaf::GetAcc(ANNpoint const& qpoint, ANNpoint& res, double /*angle2*/, ANNorthRect& bb) const
{
	double maxbox = annDist(3, bb.lo, bb.hi);
	double dist_toq = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:dist_toq)
#endif
	for (int i = 0; i < 3; ++i)
		dist_toq += (qpoint[i] - CM[i]) * (qpoint[i] - CM[i]);
	if (dist_toq < maxbox * 1e-6) //prevent self force
		return;
	double r3 = 1.0 / (dist_toq * fastsqrt(dist_toq));
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (int i = 0; i < 3; ++i)
		res[i] -= mass * (qpoint[i] - CM[i]) * r3;
	double sumQ = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:sumQ)
#endif
	for (int i = 0; i < 6; ++i)
		sumQ += std::fabs(Q[i]);
	if (sumQ < mass * dist_toq * 1e-6)
		return;
	double Qfactor = r3 / dist_toq;
	double dx = qpoint[0] - CM[0];
	double dy = qpoint[1] - CM[1];
	double dz = qpoint[2] - CM[2];
	res[0] += Qfactor * (dx * Q[0] + dy * Q[1] + dz * Q[2]);
	res[1] += Qfactor * (dx * Q[1] + dy * Q[3] + dz * Q[4]);
	res[2] += Qfactor * (dx * Q[2] + dy * Q[4] + dz * Q[5]);
	double mrr = dx * dx * Q[0] + dy * dy * Q[3] + dz * dz * Q[5] + 2 * dx * dy * Q[1] + 2 * dx * dz * Q[2] + 2 * dy * dz * Q[4];
	Qfactor *= -5 * mrr / (2 * dist_toq);
	res[0] += Qfactor * dx;
	res[1] += Qfactor * dy;
	res[2] += Qfactor * dz;
}

void ANNkd_leaf::GetAcc(std::vector<ANNpoint>& qpoint, std::vector<ANNpoint>& res, double /*angle2*/, ANNorthRect& bb) const
{
	double maxbox = annDist(3, bb.lo, bb.hi);
	size_t N = qpoint.size();
	for (size_t k = 0; k < N; ++k)
	{
		double dist_toq = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:dist_toq)
#endif
		for (int i = 0; i < 3; ++i)
			dist_toq += (qpoint[k][i] - CM[i]) * (qpoint[k][i] - CM[i]);
		if (dist_toq < maxbox * 1e-6) //prevent self force
			continue;
		double r3 = 1.0 / (dist_toq * fastsqrt(dist_toq));
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
		for (int i = 0; i < 3; ++i)
			res[k][i] -= mass * (qpoint[k][i] - CM[i]) * r3;
		double sumQ = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:sumQ)
#endif
		for (size_t i = 0; i < 6; ++i)
			sumQ += std::abs(Q[i]);
		if (sumQ < mass * dist_toq * 1e-6)
			continue;
		double Qfactor = r3 / dist_toq;
		double dx = qpoint[k][0] - CM[0];
		double dy = qpoint[k][1] - CM[1];
		double dz = qpoint[k][2] - CM[2];
		res[k][0] += Qfactor * (dx * Q[0] + dy * Q[1] + dz * Q[2]);
		res[k][1] += Qfactor * (dx * Q[1] + dy * Q[3] + dz * Q[4]);
		res[k][2] += Qfactor * (dx * Q[2] + dy * Q[4] + dz * Q[5]);
		double mrr = dx * dx * Q[0] + dy * dy * Q[3] + dz * dz * Q[5] + 2 * dx * dy * Q[1] + 2 * dx * dz * Q[2] + 2 * dy * dz * Q[4];
		Qfactor *= -5 * mrr / (2 * dist_toq);
		res[k][0] += Qfactor * dx;
		res[k][1] += Qfactor * dy;
		res[k][2] += Qfactor * dz;
	}
}

void ANNkd_leaf::GetAcc(std::vector<ANNpoint>& qpoint, std::vector<ANNpoint>& res, double /*angle2*/, ANNorthRect& bb,
	std::array<double, 4> const& /*qCM*/) const
{
	double maxbox = annDist(3, bb.lo, bb.hi);
	size_t N = qpoint.size();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (int k = 0; k < static_cast<int>(N); ++k)
	{
		double dist_toq = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:dist_toq)
#endif
		for (int i = 0; i < 3; ++i)
			dist_toq += (qpoint[k][i] - CM[i]) * (qpoint[k][i] - CM[i]);
		if (dist_toq < maxbox * 1e-6) //prevent self force
			continue;
		double r3 = 1.0 / (dist_toq * fastsqrt(dist_toq));
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
		for (int i = 0; i < 3; ++i)
			res[k][i] -= mass * (qpoint[k][i] - CM[i]) * r3;
		double sumQ = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:sumQ)
#endif
		for (int i = 0; i < 6; ++i)
			sumQ += std::abs(Q[i]);
		if (sumQ < mass * dist_toq * 1e-6)
			continue;
		double Qfactor = r3 / dist_toq;
		double dx = qpoint[k][0] - CM[0];
		double dy = qpoint[k][1] - CM[1];
		double dz = qpoint[k][2] - CM[2];
		res[k][0] += Qfactor * (dx * Q[0] + dy * Q[1] + dz * Q[2]);
		res[k][1] += Qfactor * (dx * Q[1] + dy * Q[3] + dz * Q[4]);
		res[k][2] += Qfactor * (dx * Q[2] + dy * Q[4] + dz * Q[5]);
		double mrr = dx * dx * Q[0] + dy * dy * Q[3] + dz * dz * Q[5] + 2 * dx * dy * Q[1] + 2 * dx * dz * Q[2] + 2 * dy * dz * Q[4];
		Qfactor *= -5 * mrr / (2 * dist_toq);
		res[k][0] += Qfactor * dx;
		res[k][1] += Qfactor * dy;
		res[k][2] += Qfactor * dz;
	}
}


void ANNkd_tree::GetToSend(std::vector<ANNpointArray> const& faces, std::vector<size_t>const& Nfaces, vector<ANNkd_ptr>& nodes, double angle2,
	std::vector<ANNpoint> const& normals)
{
	ANNorthRect bb(bnd_box_lo, bnd_box_hi);
	if (n_pts > 0)
		root->GetToSend(faces, Nfaces, nodes, angle2, normals, bb);
}

void ANNkd_tree::GetOpticalDepth(ANNpoint const& qpoint, std::vector<std::pair<double, double> >& res,
	double angle2) const
{
	ANNorthRect bb(bnd_box_lo, bnd_box_hi);
	root->GetOpticalDepth(qpoint, res, angle2, bb);
}

void ANNkd_tree::GetToSendOpticalDepth(std::vector<ANNpointArray> const& faces, std::vector<size_t> const& Nfaces, vector<ANNkd_ptr>& nodes, double angle2)
{
	ANNorthRect bb(bnd_box_lo, bnd_box_hi);
	std::vector<ANNorthRect> bb_faces(faces.size());
	for (size_t i = 0; i < faces.size(); ++i)
	{
		size_t N = Nfaces[i];
		bb_faces[i].hi = faces[i][0];
		bb_faces[i].lo = faces[i][0];
		for (size_t j = 1; j < N; ++j)
		{
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
			for (size_t k = 0; k < 3; ++k)
			{
				bb_faces[i].hi[k] = std::max(bb_faces[i].hi[k], faces[i][j][k]);
				bb_faces[i].lo[k] = std::min(bb_faces[i].lo[k], faces[i][j][k]);
			}
		}
	}
	if (n_pts > 0)
		root->GetToSendOpticalDepth(bb_faces, nodes, angle2, bb);
}

void ANNkd_split::GetToSend(std::vector<ANNpointArray> const& faces, std::vector<size_t>const& Nfaces, vector<ANNkd_ptr>& nodes, double angle2,
	std::vector<ANNpoint> const& normals, ANNorthRect& bb)
{
	double maxbox = annDist(3, bb.lo, bb.hi);
	double dist = DistanceToFaces(faces, Nfaces, CM, fastsqrt(maxbox / angle2), normals);
	if (dist * dist * angle2 > maxbox)
	{
		nodes.push_back(this);
		return;
	}

	double lv = bb.lo[cut_dim];
	double hv = bb.hi[cut_dim];

	if (child[1] != KD_TRIVIAL)
	{
		bb.lo[cut_dim] = cut_val;
		child[1]->GetToSend(faces, Nfaces, nodes, angle2, normals, bb);
		bb.lo[cut_dim] = lv;
	}
	if (child[0] != KD_TRIVIAL)
	{
		bb.hi[cut_dim] = cut_val;
		child[0]->GetToSend(faces, Nfaces, nodes, angle2, normals, bb);
		bb.hi[cut_dim] = hv;
	}
}

void ANNkd_leaf::GetToSend(std::vector<ANNpointArray> const& /*faces*/, std::vector<size_t>const& /*Nfaces*/, vector<ANNkd_ptr>& nodes, double /*angle2*/,
	std::vector<ANNpoint> const& /*normals*/, ANNorthRect&/*bb*/)
{
	nodes.push_back(this);
}
