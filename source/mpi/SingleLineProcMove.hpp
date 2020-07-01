/*! \file SingleLineProcMove.hpp
\brief A class that tries to maintain a constant number of points by placing cpus in a line on the x axis
\author Elad Steinberg
*/

#ifndef SINGLELINEPROC_HPP
#define SINGLELINEPROC_HPP 1
#include "ProcessorUpdate3D.hpp"

//! \brief A load balancing scheme aiming for the same number of points in each process
class SingleLineProcMove : public ProcessorUpdate3D
{
public:
	void Update(Tessellation3D& tproc, Tessellation3D const& tlocal)const;
};
#endif //SINGLELINEPROC_HPP
