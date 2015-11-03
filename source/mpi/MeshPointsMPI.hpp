/*! \file MeshPointsMPI.hpp
  \brief Generates grids for parallel simulations
  \author Elad Steinberg
 */

#ifndef MESHPOINTSMPI_HPP
#define MESHPOINTSMPI_HPP 1
#ifdef RICH_MPI
#define _USE_MATH_DEFINES
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "../tessellation/ConvexHull.hpp"
#include "../tessellation/VoronoiMesh.hpp"
#include "../misc/utils.hpp"

/*!
\brief Creates a cartesian mesh for MPI
\param nx The total number of point in the x direction
\param ny The total number of point in the y direction
\param tess The tessellation of the processors
\param lower_left Lower left point
\param upper_right Upper right point
\return The set of points corresponding to the local process
*/
vector<Vector2D> SquareMeshM(int nx,int ny,Tessellation const& tess,
			     Vector2D const& lower_left,
			     Vector2D const& upper_right);

/*!
  \brief Generates a round grid with constant point density
  \param PointNum The number of points.
  \param Rmin The min radius
  \param Rmax The max radius
  \param bottomleft The lower left corner of a boundaing box to cut off the circle
  \param topright The top right corner of a boundaing box to cut off the circle
  \param tess The tessellation of the processors
  \param xc X of circle center
  \param yc Y of circle center
  \return List of two dimensional points
*/
vector<Vector2D> CirclePointsRmaxM(int PointNum,double Rmin,double Rmax,
	Vector2D const& bottomleft,Vector2D const& topright,
	Tessellation const& tess,double xc=0,double yc=0);
/*!
  \brief Generates a round grid with r^alpha point density
  \param PointNum The number of points.
  \param Rmin The min radius
  \param Rmax The max radius
  \param xc X of circle center
  \param yc Y of circle center
  \param alpha The point density, should not be -1 or -2
  \param tess The tessellation of the processors
  \return List of two dimensional points
*/
vector<Vector2D> CirclePointsRmax_aM(int PointNum,double Rmin,double Rmax,
	double xc,double yc,double alpha,Tessellation const& tess);

/*!
  \brief Creates a circle of evenly spaced points
  \param point_number Number of points along the circumference
  \param radius Radius of the circle
  \param center Position of the center of the circle
  \param tproc The tessellation of the processors
  \return List of two dimensional points
*/
vector<Vector2D> circle_circumferenceM(int point_number,double radius,
	Vector2D const& center,Tessellation const& tproc);

/*!
\brief Creates a uniform random mesh for MPI
\param npoints The total number of point
\param tess The tessellation of the processors
\param lowerleft The lower left point of the domain
\param upperright The upper right point of the domain
\return The set of points corresponding to the local process
*/
vector<Vector2D> RandSquare(int npoints,Tessellation const& tess,
	Vector2D const& lowerleft,Vector2D const& upperright);

#endif
#endif //MESHPOINTSMPI_HPP
