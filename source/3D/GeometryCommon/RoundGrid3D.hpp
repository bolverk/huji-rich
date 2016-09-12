/*! \file RoundGrid3D.hpp
  \brief Makes the initial cells rounder
  \author Elad Steinberg
 */

#ifndef ROUNDGRID3D
#define ROUNDGRID3D 1

#include "Voronoi3D.hpp"

/*!
	\brief Makes the cells rounder
	\param points The initial points
	\param ll The lower left corner of the domain
	\param ur The upper right corner of the domain
	\param NumberIt The number of correction iterations
	\param tess The tessellation
*/
#ifdef RICH_MPI
//!	\param tproc The tessellation of processors
#endif
/*!
	\return The points that give a rounder tessellation
*/
vector<Vector3D> RoundGrid3D(vector<Vector3D> const& points,Vector3D const& ll,Vector3D const& ur,
	size_t NumberIt=10,
	#ifdef RICH_MPI
		   Tessellation3D const* tproc=0,
	#endif
	Tessellation3D *tess=0);
#endif //ROUNDGRID3D
