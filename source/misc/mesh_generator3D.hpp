/*! \file mesh_generator3D.hpp
\brief Set of functions to generate points.
\author Elad Steinberg
*/
#ifndef MESHGENERATOR3D_HPP
#define MESHGENERATOR3D_HPP 1

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif // _MSC_VER
#include <vector>
#include <cmath>
#include "../3D/GeometryCommon/Voronoi3D.hpp"
#include <algorithm>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

/*! \brief Generates a cartesian mesh
\param nx Number of points along the x axis
\param ny Number of points along the y axis
\param nz Number of points along the z axis
\param lower_left Lower left point
\param upper_right Upper right point
\return Set of three dimensional points
*/
vector<Vector3D> CartesianMesh(std::size_t nx, std::size_t ny, std::size_t nz, Vector3D const& lower_left,
	Vector3D const& upper_right);

/*!
\brief Generates a random grid with uniform point density and a constant seed
\param PointNum The number of points.
\param ll The lower left point of the domain
\param ur The upper right point of the domain
\return List of three dimensional points
*/
vector<Vector3D> RandRectangular(std::size_t PointNum, Vector3D const& ll, Vector3D const& ur,Voronoi3D const* tproc = 0);

vector<Vector3D> RandRectangular(std::size_t PointNum, Vector3D const& ll, Vector3D const& ur, boost::mt19937_64 &gen);

vector<Vector3D> RandSphereR(std::size_t PointNum, Vector3D const& ll, Vector3D const& ur, double Rmin, double Rmax,
	Vector3D center = Vector3D(),Voronoi3D const* tproc = 0);

vector<Vector3D> RandSphereR2(std::size_t PointNum, Vector3D const& ll, Vector3D const& ur,double Rmin,double Rmax
	, Vector3D center = Vector3D(), Voronoi3D const* tproc = 0);

vector<Vector3D> RandSphereR1(std::size_t PointNum, Vector3D const& ll, Vector3D const& ur, double Rmin, double Rmax,
	Vector3D center = Vector3D(),Voronoi3D const* tproc = 0);

vector<Vector3D> RandSphereRa(std::size_t PointNum, Vector3D const& ll, Vector3D const& ur, double Rmin, double Rmax,double a, Vector3D const& center,
	Voronoi3D const* tproc = 0);

#ifdef RICH_MPI
/*!
\brief Generates a random grid with uniform point density and a constant seed
\param PointNum The total number of points to be in all cpus combined.
\param tproc The tessellation of the processors
\return List of three dimensional points
*/
vector<Vector3D> RandPointsMPI(Voronoi3D const& tproc, size_t PointNum);
#endif

#endif //MESHGENERATOR3D_HPP

