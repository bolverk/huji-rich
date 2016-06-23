/*! \file mesh_generator.hpp
  \brief Set of functions to generate points.
  \author Elad Steinberg
*/
#ifndef MESHGENERATOR_HPP
#define MESHGENERATOR_HPP 1

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif // _MSC_VER
#include <vector>
#include <cmath>
#include "../tessellation/geometry.hpp"
#include <algorithm>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

/*! \brief Generates a cartesian mesh
  \param nx Number of points along the x axis
  \param ny Number of points along the y axis
  \param lower_left Lower left point
  \param upper_right Upper right point
  \return Set of two dimensional points
*/
std::vector<Vector2D> cartesian_mesh(int nx, int ny,
				Vector2D const& lower_left,
				Vector2D const& upper_right);

/*!
  \brief Generates a round grid with constant point density
  \param PointNum The number of points.
  \param Rmin The min radius
  \param Rmax The max radius
  \param xc X of circle center
  \param yc Y of circle center
  \param xmin Left edge of confining rectangle
  \param xmax Right edge of confining rectangle
  \param ymax Upper edge of confining rectangle
  \param ymin Lower edge of confining rectangle
  \return List of two dimensional points
*/
std::vector<Vector2D> CirclePointsRmax(int PointNum,double Rmin,double Rmax,
	double xmax, double ymax, double xmin, double ymin, double xc = 0, double yc = 0);

/*!
  \brief Generates a round grid with 1/r^2 point density confined to a rectangle given by xmin,xmax,ymin and ymax.
  \param PointNum The number of points.
  \param Rmin The min radius
  \param Rmax The max radius
  \param xc X of circle center
  \param yc Y of circle center
  \param xmin Left edge of confining rectangle
  \param xmax Right edge of confining rectangle
  \param ymax Upper edge of confining rectangle
  \param ymin Lower edge of confining rectangle
  \return List of two dimensional points
*/
std::vector<Vector2D> CirclePointsRmax_2(int PointNum,double Rmin,double Rmax,
				    double xc=0,double yc=0,double xmax=1,double ymax=1,double xmin=-1,
				    double ymin=1);
/*!
  \brief Generates a round grid with 1/r point density confined to a rectangle given by xmin,xmax,ymin and ymax.
  \param PointNum The number of points.
  \param Rmin The min radius
  \param Rmax The max radius
  \param xc X of circle center
  \param yc Y of circle center
  \param xmin Left edge of confining rectangle
  \param xmax Right edge of confining rectangle
  \param ymax Upper edge of confining rectangle
  \param ymin Lower edge of confining rectangle
  \return List of two dimensional points
*/
std::vector<Vector2D> CirclePointsRmax_1(int PointNum,double Rmin,double Rmax,
				    double xc=0,double yc=0,double xmax=1,double ymax=1,double xmin=-1,double ymin=-1);

/*!
  \brief Creates a circle of evenly spaced points
  \param point_number Number of points along the circumference
  \param radius Radius of the circle
  \param center Position of the center of the circle
  \return List of two dimensional points
*/
std::vector<Vector2D> circle_circumference(size_t point_number,
					   double radius,
					   Vector2D const& center);

/*!
  \brief Creates a line of evenly spaced points y=slope*x+b
  \param PointNum The number of points
  \param xmin The minimum x of the line
  \param xmax The maximum x of the line
  \param ymin The minimum y of the line
  \param ymax The maximum y of the line
  \return List of two dimensional points
*/
std::vector<Vector2D> Line(int PointNum,double xmin,double xmax,double ymin,double ymax);

/*!
  \brief Generates a round grid with r^alpha point density confined to a rectangle given by xmin,xmax,ymin and ymax.
  \param PointNum The number of points.
  \param Rmin The min radius
  \param Rmax The max radius
  \param xc X of circle center
  \param yc Y of circle center
  \param xmin Left edge of confining rectangle
  \param xmax Right edge of confining rectangle
  \param ymax Upper edge of confining rectangle
  \param ymin Lower edge of confining rectangle
  \param alpha The point density, should not be -1 or -2
  \return List of two dimensional points
*/
std::vector<Vector2D> CirclePointsRmax_a(int PointNum,double Rmin,double Rmax,
				    double xc,double yc,double xmax,double ymax,double xmin,double ymin,
				    double alpha);

/*!
  \brief Generates a rectangular grid with random 1/r point density
  \param PointNum The number of points.
  \param xl The left boundary
  \param xr The right boundary
  \param yd The lower boundary
  \param yu The upper boundary
  \param minR The inner radius in which there are no points
  \param xc The X of center of the circle
  \param yc The Y of center of the circle
  \return List of two dimensional points
*/
std::vector<Vector2D> RandPointsR(int PointNum,double xl=-0.5,double xr=0.5,
			     double yd=-0.5,double yu=0.5,double minR=0,double xc=0,
				 double yc=0);

/*!
  \brief Generates a random rectangular grid with uniform point density and a constant seed
  \param PointNum The number of points.
  \param lowerleft The lower left point of the domain
  \param upperright The upper right point of the domain
  \return List of two dimensional points
*/
std::vector<Vector2D> RandSquare(int PointNum,Vector2D const& lowerleft, Vector2D const& upperright);


/*!
\brief Generates a random rectangular grid with uniform point density and a constant seed
\param PointNum The number of points.
\param xl The left boundary
\param xr The right boundary
\param yd The lower boundary
\param yu The upper boundary
\return List of two dimensional points
*/
std::vector<Vector2D> RandSquare(int PointNum, double xl = -0.5, double xr = 0.5,
	double yd = -0.5, double yu = 0.5);

/*!
  \brief Generates a random rectangular grid with uniform point density. This is when reseting the seed between calls isn't wanted
  \param PointNum The number of points.
  \param eng The random number generator
  \param xl The left boundary
  \param xr The right boundary
  \param yd The lower boundary
  \param yu The upper boundary
  \return List of two dimensional points
*/

std::vector<Vector2D> RandSquare(int PointNum,boost::random::mt19937 &eng,
	double xl=-0.5,double xr=0.5,double yd=-0.5,double yu=0.5);

/*!
  \brief Generates a random round grid with r^(1-a) point density
  \param PointNum The number of points.
  \param Rmin The min radius
  \param Rmax The max radius
  \param alpha The radial density of the points, shouldn't be 1
  \param lowerleft The lowerleft corner of the domain
  \param upperright The upperright corner of the domain
  \param center The center of the coordinates system
  \return List of two dimensional points
*/

std::vector<Vector2D> RandPointsRa(int PointNum,double Rmin,double Rmax,double alpha,
	Vector2D const& lowerleft,Vector2D const& upperright,Vector2D const& center=Vector2D(0,0));

/*!
  \brief Generates a random round grid with 1/r point density
  \param PointNum The number of points.
  \param Rmin The min radius
  \param Rmax The max radius
  \param xc X of circle center
  \param yc Y of circle center
  \return List of two dimensional points
*/
std::vector<Vector2D> RandPointsRmax(int PointNum,double Rmin,double Rmax,
				double xc=0,double yc=0);

/*!
\brief Generates a random round grid with 1/y point density where y is the cylindrical r coordinate
\param PointNum The number of points.
\param ll The lower left point of the domain. Should be ll.y>0
\param ur The upper right point of the domain
\return List of two dimensional points
*/
std::vector<Vector2D> RandPointsCylinder(int PointNum, Vector2D const& ll,Vector2D const& ur);

#endif //MESHGENERATOR_HPP
