/*! \file RoundGrid.hpp
  \brief Makes the initial cells rounder
  \author Elad Steinberg
 */

#ifndef ROUNDGRID
#define ROUNDGRID 1
#include "VoronoiMesh.hpp"

#define DEFAULT_ITER_NUM 10

#ifdef RICH_MPI
vector<Vector2D> RoundGridV
(vector<Vector2D> const& points,
 const OuterBoundary& bc,
 int NumberIt=DEFAULT_ITER_NUM);
#endif // RICH_MPI

vector<Vector2D> RoundGridV
(vector<Vector2D> const& points,
 const OuterBoundary& bc,
#ifdef RICH_MPI
 const Tessellation& meta,
#endif
 int NumberIt=DEFAULT_ITER_NUM);

/*!
	\brief Makes the cells rounder
	\param points The initial points
	\param bc The outer boundary conditions
	\param NumberIt The number of correction iterations
	\param tess The tessellation
*/
#ifdef RICH_MPI
//!	\param tproc The tessellation of processors
#endif
/*!
	\return The points that give a rounder tessellation
*/
vector<Vector2D> RoundGrid
(vector<Vector2D> const& points,
 const OuterBoundary& bc,
 Tessellation& tess,
#ifdef RICH_MPI
 const Tessellation& meta,
#endif
 int NumberIt=DEFAULT_ITER_NUM);

#endif //ROUNDGRID
