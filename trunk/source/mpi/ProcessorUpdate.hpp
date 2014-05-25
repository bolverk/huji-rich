/*!
  \brief Abstract class for motion of the processor points
  \author Elad Steinberg
*/
#ifndef PROCUPDATE
#define PROCUPDATE 1

#define _USE_MATH_DEFINES 
#include <cmath>
#include "../newtonian/two_dimensional/OuterBoundary.hpp"
#include "../tessellation/tessellation.hpp"
#include "mpi_macro.hpp"

//! \brief Updates the positions of the processes
class ProcessorUpdate
{
public:
  /*!
    \brief Moves the processor tessellation
    \param tproc The tessellation of the processors
    \param tlocal The tessellation of the local mesh points
  */
  virtual void Update(Tessellation &tproc,Tessellation const& tlocal)const=0;

  //! \brief virtual destructor
  virtual ~ProcessorUpdate(void);
};


#endif //PROCUPDATE
