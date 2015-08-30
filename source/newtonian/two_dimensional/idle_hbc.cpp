#include "idle_hbc.hpp"

IdleHBC::IdleHBC(void) {}

vector<pair<size_t,Extensive> > IdleHBC::operator()
  (const Tessellation& /*tessellation*/,
   const vector<ComputationalCell>& /*cells*/) const
{
  return vector<pair<size_t,Extensive> >();
}
