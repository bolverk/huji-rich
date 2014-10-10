/*! \file RoundGrid.hpp
  \brief Makes the initial cells rounder
  \author Elad Steinberg
 */

#ifndef ROUNDGRID
#define ROUNDGRID 1
#include "VoronoiMesh.hpp"

/*!
	\brief Makes the cells rounder
	\param points The initial points
	\param bc The outer boundary conditions
	\param NumberIt The number of correction iterations
	\param InnerNum The number of first points not to move
	\param tess The tessellation
	\param tproc The tessellation of processors
	\return The points that give a rounder tessellation
*/
vector<Vector2D> RoundGrid(vector<Vector2D> const& points,
	OuterBoundary const* bc,int NumberIt=10,
	int InnerNum=0,
			   #ifdef RICH_MPI
			   Tessellation const* tproc=0,
			   #endif
			   Tessellation *tess=0);
#endif //ROUNDGRID
