/*! \file ProcessorUpdate.hpp
\brief Abstract class for motion of the processor points
\author Elad Steinberg
*/
#ifndef PROCUPDATE
#define PROCUPDATE 1

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif // _MSC_VER
#include <cmath>
#include <algorithm>
#include "../newtonian/two_dimensional/OuterBoundary.hpp"
#include "../tessellation/tessellation.hpp"

//! \brief Updates the positions of the processes
class ProcessorUpdate
{
public:
#ifdef RICH_MPI
	/*!
	\brief Moves the processor tessellation
	\param tproc The tessellation of the processors
	\param tlocal The tessellation of the local mesh points
	*/
	virtual void Update(Tessellation &tproc,Tessellation const& tlocal)const=0;
#endif
	//! \brief virtual destructor
	virtual ~ProcessorUpdate(void);
#ifdef RICH_MPI
	/*!
	\brief Calcualtes the load imbalance as max(number of points per proc)/(avg per proc)
	\param tlocal The local tesselaltion
	\param total The total number of points, given as output
	\return The load imbalance
	*/
	double GetLoadImbalance(Tessellation const& tlocal,int &total)const;
#endif
};

#endif //PROCUPDATE
