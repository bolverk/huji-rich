#ifndef SCALAR_INTERPOLATION_HPP
#define SCALAR_INTERPOLATION_HPP 1

#include "../../tessellation/tessellation.hpp"

enum SInterpolationType{SInBulk,SBoundary};

//! \brief Abstract class for interpolating scalars
class ScalarInterpolation
{
public:
	/*! \brief Prepares the interpolation by calculating all of the slopes
	\param tessellation The tessellation
	\param tracers The passive tracers
	\param dt The time step
	\param time The simulation time
	*/
	virtual void Prepare(Tessellation const* tessellation,
		vector<vector<double> > const& tracers,
		double dt,double time) = 0;
	/*!
	\brief Interpolates the passive scalars
	\param tessellation The tessellation
	\param tracers The passive tracers
	\param dt The time step
	\param edge The edge to do the interpolation on
	\param side The edge side to interpolate
	\param interptype Interpolation type (in the bulk or boundary)
	\return The interpolated scalars
	*/
	virtual vector<double> Interpolate(Tessellation const* tessellation,
		vector<vector<double> > const& tracers,double dt,
		Edge const& edge,int side,SInterpolationType interptype)
		const = 0;
	/*!
	\brief Virtual destructor
	*/
  virtual ~ScalarInterpolation(void);
};

#endif // SCALAR_INTERPOLATION_HPP
