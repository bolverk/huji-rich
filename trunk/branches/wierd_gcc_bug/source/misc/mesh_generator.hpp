/*!
	\brief Set of functions to generate points.
	\author Elad Steinberg
*/
#ifndef MESHGENERATOR_HPP
#define MESHGENERATOR_HPP 1
#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include "../tessellation/geometry.hpp"
#include <algorithm>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
typedef boost::mt19937_64 base_generator_type;

using namespace std;

/*!
	\brief Generates a square grid around (0,0)
	\param sidex The x length
	\param sidey The y length
	\param nx The number of points in the x direction
	\param ny The number of points in the y direction
	\param centerd If true then the mesh is centered around (0,0) else (0,0) is the lower left point
*/

vector<Vector2D> SquareMesh(int nx,int ny,double sidex=1,double sidey=1,
	bool centerd=true);

/*! \brief Generates a cartesian mesh
  \param nx Number of points along the x axis
  \param ny Number of points along the y axis
  \param lower_left Lower left point
  \param upper_right Upper right point
 */
vector<Vector2D> cartesian_mesh(int nx, int ny,
				Vector2D const& lower_left,
				Vector2D const& upper_right);


/*!
	\brief Generates a round grid with constant point density
	\param PointNum The number of points.
	\param Rmin The min radius
	\param Rmax The max radius
	\param xc X of circle center
	\param yc Y of circle center
*/

vector<Vector2D> CirclePointsRmax(int PointNum,double Rmin,double Rmax,
	double xc=0,double yc=0);

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
*/
vector<Vector2D> CirclePointsRmax_2(int PointNum,double Rmin,double Rmax,
	double xc=0,double yc=0,double xmax=1,double ymax=1,double xmin=0,
	double ymin=0);
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
*/
vector<Vector2D> CirclePointsRmax_1(int PointNum,double Rmin,double Rmax,
	double xc=0,double yc=0,double xmax=1,double ymax=1,double xmin=0,double ymin=0);

/*!
	\brief Creates a circle of evenly spaced points
	\param PointNum The number of points
	\param R The radius of the circle
	\param xc The X center
	\param yc The Y center
*/
vector<Vector2D> Circle(int PointNum,double R,double xc=0,double yc=0);

/*!
	\brief Creates a line of evenly spaced points y=slope*x+b
	\param PointNum The number of points
	\param xmin The minimum x of the line
	\param xmax The maximum x of the line
	\param ymin The minimum y of the line
	\param ymax The maximum y of the line
*/
vector<Vector2D> Line(int PointNum,double xmin,double xmax,double ymin,double ymax);

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
*/

vector<Vector2D> CirclePointsRmax_a(int PointNum,double Rmin,double Rmax,
	double xc,double yc,double xmax,double ymax,double xmin,double ymin,
	double alpha);

/*!
	\brief Generates a square grid around 0,0 with small pertubations
	\param sidex The x length
	\param sidey The y length
	\param nx The number of points in the x direction
	\param ny The number of points in the y direction
	\param mag The magnitude of the small pertubations, should be less than 0.5
*/

vector<Vector2D> SquarePertubed(int nx,int ny,double sidex=1,double sidey=1,
	double mag=0.01);

/*!
	\brief Generates a square grid with 1/r point density
	\param PointNum The number of points.
	\param xl The left boundary
	\param xr The right boundary
	\param yd The lower boundary
	\param yu The upper boundary
	\param minR The inner radius in which there are no points
*/

vector<Vector2D> RandPointsR(int PointNum,double xl=-0.5,double xr=0.5,
	double yd=-0.5,double yu=0.5,double minR=0);
/*!
	\brief Generates a square grid with unifrom point density
	\param PointNum The number of points.
	\param xl The left boundary
	\param xr The right boundary
	\param yd The lower boundary
	\param yu The upper boundary
*/

vector<Vector2D> RandSquare(int PointNum,double xl=-0.5,double xr=0.5,
	double yd=-0.5,double yu=0.5);
/*!
	\brief Generates a round grid with 1/r point density
	\param PointNum The number of points.
	\param Rmin The min radius
	\param Rmax The max radius
	\param xc X of circle center
	\param yc Y of circle center
*/

vector<Vector2D> RandPointsRmax(int PointNum,double Rmin,double Rmax,
	double xc=0,double yc=0);


#endif //MESHGENERATOR_HPP
