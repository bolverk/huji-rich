/*! \file trigger.hpp
  \brief Trigger for diagnostic function
  \author Almog Yalinewich
 */

#ifndef TRIGGER_HPP
#define TRIGGER_HPP 1

#include "../two_dimensional/hdsim2d.hpp"

class Trigger
{
public:

  virtual bool operator()(const hdsim& sim) = 0;

  virtual ~Trigger(void);
};

class ConstantTimeInterval: public Trigger
{
public:
  ConstantTimeInterval(double dt, double t_next=0);

  bool operator()(const hdsim& sim);

private:

  const double dt_;
  mutable double t_next_;
};

#endif // TRIGGER_HPP
