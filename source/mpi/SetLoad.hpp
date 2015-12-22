/*! \file SetLoad.hpp
  \brief Function for setting the load balance based on equal cells per processor
  \author Elad Steinberg
*/
#ifndef SETLOAD
#define SETLOAD 1

#include "ConstNumberPerProc.hpp"
#include "../tessellation/VoronoiMesh.hpp"
#include "../newtonian/two_dimensional/diagnostics.hpp"
#include "../misc/simple_io.hpp"

/*!
	\brief Corrects the load between processors based on number of cells per processor
	\param tproc The tessellation of the processors
	\param points The local cell points of the current rank (the points for this current cpu), may be redistributed among other cpus
	\param outer The geometric boundary conditions
	\param Niter The number of correction iterations to use
	\param speed How fast to make the correction each iteration in units of cpu cell size
	\param mode The operating mode, 1=Hybrid, 2=Pressure based, 3= Density based
*/
void SetLoad(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
	outer,int Niter=100,double speed=0.04,int mode=2);

/*!
	\brief Corrects the load between processors based on number of cells per processor, continues until target load is reached
	\param tproc The tessellation of the processors
	\param points The local cell points of the current rank (the points for this current cpu), may be redistributed among other cpus
	\param outer The geometric boundary conditions
	\param TargetLoad The target load balance, must be higher than 1, recommened to be 1.5
	\param speed How fast to make the correction each iteration in units of cpu cell size
	\param mode The operating mode, 1=Hybrid, 2=Pressure based, 3= Density based
*/
void SetLoad(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
	outer,double TargetLoad,double speed=0.04,int mode=2);

#endif //SETLOAD
