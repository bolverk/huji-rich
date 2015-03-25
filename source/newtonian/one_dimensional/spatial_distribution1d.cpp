/*!  \brief Abstract class for initial conditions
  \author Almog Yalinewich
 */

#include "spatial_distribution1d.hpp"

SpatialDistribution1D::~SpatialDistribution1D(void) {}

Uniform::Uniform(double val): _val(val) {}

double Uniform::operator()(double /*x*/) const
{
  return _val;
}

Step::Step(double val1, double val2,
	   double steppos):
  _val1(val1), _val2(val2), _steppos(steppos) {}

double Step::operator()(double x) const
{
  if (x>_steppos)
    return _val2;
  else
    return _val1;
}

TwoSteps::TwoSteps(double ivl, double ilip,
		   double ivm, double irip,
		   double ivr):
  vl(ivl),
  vr(ivr),
  vm(ivm),
  lip(ilip),
  rip(irip) {}

double TwoSteps::operator()(double x) const
{
  if (x<lip)
    return vl;
  else if (x>rip)
    return vr;
  else
    return vm;
}
