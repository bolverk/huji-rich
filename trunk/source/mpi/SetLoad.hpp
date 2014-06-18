/*!
  \brief Function for setting the load balance based on equal cells per processor
  \author Elad Steinberg
*/
#ifndef SETLOAD
#define SETLOAD 1

#include "ConstNumberPerProc.hpp"
#include "../tessellation/VoronoiMesh.hpp"
#include "../newtonian/two_dimensional/diagnostics.hpp"

/*!
	\brief Corrects the load between processors based on number of cells per processor
	\param tproc The tessellation of the processors
	\param points The local cell points of the current rank (the points for this current cpu)
	\param outer The geometric boundary conditions
	\param Nbest The ideal number of points per processor
	\param Niter The number of correction iterations to use
*/
void SetLoad(Tessellation &tproc,vector<Vector2D> const& points,OuterBoundary const&
	outer,int Nbest,int Niter=10);

#endif //SETLOAD