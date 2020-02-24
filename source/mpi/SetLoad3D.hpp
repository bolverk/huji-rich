/*! \file SetLoad3D.hpp
\brief Function for setting the load balance based on equal cells per processor
\author Elad Steinberg
*/
#ifndef SETLOAD3D
#define SETLOAD3D 1

#ifdef RICH_MPI

#include "mpi_commands.hpp"
#include "ConstNumberPerProc3D.hpp"
#include "../newtonian/three_dimensional/computational_cell.hpp"
#include "../3D/GeometryCommon/Voronoi3D.hpp"

/*!
\brief Corrects the load between processors based on number of cells per processor
\param tproc The tessellation of the processors
\param points The local mesh points of the current rank (the points for this current cpu), may be redistributed among other cpus
\param Niter The number of correction iterations to use
\param speed How fast to make the correction each iteration in units of cpu cell size
\param mode The operating mode, 1=Hybrid, 2=Pressure based, 3= Density based
\param round The factor to enhance the cells rounding mechanisim
*/
void SetLoad(Voronoi3D &tproc, vector<Vector3D> &points,size_t Niter = 300, double speed = 0.03, int mode = 2,
	double round = 0.05, bool display = false);
/*!
\brief Corrects the load between processors based on number of cells per processor
\param tproc The tessellation of the processors
\param cells he local computational cell points of the current rank (the points for this current cpu), may be redistributed among other cpus
\param points The local mesh points of the current rank (the points for this current cpu), may be redistributed among other cpus
\param Niter The number of correction iterations to use
\param speed How fast to make the correction each iteration in units of cpu cell size
\param mode The operating mode, 1=Hybrid, 2=Pressure based, 3= Density based
\param round The factor to enhance the cells rounding mechanisim
*/

void SetLoad(Voronoi3D &tproc, vector<Vector3D> &points,vector<ComputationalCell3D> &cells, size_t Niter = 300, double speed = 0.03, int mode = 21,
	double round = 0.05, bool display = false);


#endif // RICH_MPI

#endif //SETLOAD3D