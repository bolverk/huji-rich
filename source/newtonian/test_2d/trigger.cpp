#include "trigger.hpp"

ConstantTimeInterval::ConstantTimeInterval(double dt, double t_next):
  dt_(dt), t_next_(t_next) {}

bool ConstantTimeInterval::operator()(const hdsim& sim)
{
  if(sim.GetTime()>t_next_){
    t_next_ += dt_;
    return true;
  }
  else
    return false;
}
