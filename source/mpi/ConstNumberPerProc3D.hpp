/*! \file ConstNumberPerProc3D.hpp
\brief A class that tries to maintain a constant number of points per processor by solving eq 68 in AREPO's paper
\author Elad Steinberg
*/

#ifndef CONSTPERPROC3D
#define CONSTPERPROC3D 1
#include "ProcessorUpdate3D.hpp"

//! \brief A load balancing scheme aiming for the same number of points in each process
class ConstNumberPerProc3D : public ProcessorUpdate3D
{
public:
	/*!
	\brief Class constructor
	\param speed The maximum change (in cell radii) of the processor movement
	\param RoundSpeed By which factor is the rounding mechanisim larger than the movement speed of the mesh points
	\param mode The operating mode, 1=Hybrid, 2=Pressure based, 3= Density based
	\param Hilbert Flag to use hilbert ordering for load balance
	*/
	ConstNumberPerProc3D(double speed = 0.01,double RoundSpeed = 0.025, int mode = 2,bool Hilbert=false);

	/*!
	\brief Updates the load balance, does one iteration
	\param tproc The processors tessellation.
	\param tlocal The local tessellation
	*/
	void Update(Tessellation3D &tproc, Tessellation3D const& tlocal)const;

	//! \brief Class destructor
	~ConstNumberPerProc3D(void);

private:
	const double speed_;
	const double RoundSpeed_;
	const int mode_;
	const bool Hilbert_;
	mutable int run_counter_;
};
#endif //CONSTPERPROC3D