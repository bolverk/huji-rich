/*! \file ConstNumberPerProc.hpp
  \brief A class that tries to maintain a constant number of points per processor by solving eq 68 in AREPO's paper
  \author Elad Steinberg
 */

#ifndef CONSTPERPROC
#define CONSTPERPROC 1
#include "ProcessorUpdate.hpp"

//! \brief A load balancing scheme aiming for the same number of points in each process
class ConstNumberPerProc: public ProcessorUpdate
{
public:
	/*!
	\brief Class constructor
	\param outer The outer geometric boundary conditions
	\param speed The maximum change (in cell radii) of the processor movement
	\param RoundSpeed By which factor is the rounding mechanisim larger than the movement speed of the mesh points
	\param mode The operating mode, 1=Hybrid, 2=Pressure based, 3= Density based
	*/
	ConstNumberPerProc(OuterBoundary const& outer,double speed=0.03,
		double RoundSpeed=2,int mode=2);

	/*!
	\brief Updates the load balance, does one iteration
	\param tproc The processors tessellation.
	\param tlocal The local tessellation
	*/
	void Update(Tessellation &tproc,Tessellation const& tlocal)const;

	//! \brief Class destructor
	~ConstNumberPerProc(void);

private:
	OuterBoundary const& outer_;
	const double speed_;
	const double RoundSpeed_;
	const int mode_;
};
#endif //CONSTPERPROC
