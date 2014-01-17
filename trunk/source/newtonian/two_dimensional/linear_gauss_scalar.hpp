/*! \file linear_gauss_scalar.hpp
  \brief Linear interpolation for the tracers
  \author Elad Steinberg
  \deprecated Due to merging of hydro and tracer interpolations
 */

#ifndef LINEAR_GAUSS_SCALAR_HPP
#define LINEAR_GAUSS_SCALAR_HPP 1

#include "scalar_interpolation.hpp"
#include "spatial_reconstruction.hpp"
#include "HydroBoundaryConditions.hpp"
#include "CustomEvolution.hpp"
#include "../../misc/utils.hpp"
#include "../../tessellation/geometry.hpp"
#include "../../misc/universal_error.hpp"
//! \brief Linear Gauss interpolation for scalar fields
class LinearGaussScalar: public ScalarInterpolation
{
public:
	/*!
	\brief Class constructor
	\param sr The spatial reconstruction for the hydro cells
	\param hbc The hydro boundary conditions
	\param theta The parameter that determines the slope limiter strength
	*/
	LinearGaussScalar(SpatialReconstruction const* sr,HydroBoundaryConditions
		const* hbc,double theta=0.5);

	/*!
	\brief Class destructor
	*/
	~LinearGaussScalar(void);

	void Prepare(Tessellation const* tessellation,
		vector<vector<double> > const& tracers,
		double dt,double time);

	vector<double> Interpolate(Tessellation const* tessellation,
		vector<vector<double> > const& tracers,double dt,
		Edge const& edge,int side,SInterpolationType interptype)const;

private:
	SpatialReconstruction const* sr_;
	vector<vector<Vector2D> > slopes_;
	HydroBoundaryConditions const* hbc_;
	double theta_;

  LinearGaussScalar(const LinearGaussScalar& origin);
  LinearGaussScalar& operator=(const LinearGaussScalar& origin);
};

#endif // LINEAR_GAUSS_SCALAR_HPP
