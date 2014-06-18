#ifndef CONSTPERPROC
#define CONSTPERPROC 1
#include "ProcessorUpdate.hpp"
//! \brief A class that tries to maintain a constant number of points per processor by solving eq 68 in AREPO's paper
class ConstNumberPerProc: public ProcessorUpdate
{
public:
	/*!
	\brief Class constructor
	\param outer The outer geometric boundary conditions
	\param npercell The ideal number of points per processor
	\param speed The maximum change (in cell radii) of the processor movement
	*/
	ConstNumberPerProc(OuterBoundary const& outer,int npercell,double speed=0.05);

	void Update(Tessellation &tproc,Tessellation const& tlocal)const;

	//! \brief Class destructor
	~ConstNumberPerProc(void);
private:
	OuterBoundary const& outer_;
	const double PointsPerProc_;
	double speed_;
};
#endif //CONSTPERPROC