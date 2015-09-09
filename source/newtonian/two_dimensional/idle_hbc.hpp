#ifndef IDLE_HBC_HPP
#define IDLE_HBC_HPP 1

#include "HydroBoundaryConditions.hpp"

//! \brief Does not compute boundary conditions at all (return an empty list)
class IdleHBC: public HydroBoundaryConditions
{
public:

  IdleHBC(void);

  vector<pair<size_t,Extensive> > operator()
  (const Tessellation& tessellation,
   const vector<ComputationalCell>& cells) const;
};

#endif // IDLE_HBC_HPP
