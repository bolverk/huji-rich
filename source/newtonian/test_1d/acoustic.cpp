//! \file acoustic.cpp

#include <cmath>
#include "acoustic.hpp"
#include "../../misc/universal_error.hpp"

AcousticInitCond::AcousticInitCond(double d0, double p0,
				   EquationOfState const& eos,
				   double dd, double wavelength):
  c0_(eos.dp2c(d0,p0)),
  density_(dd,wavelength,0,d0),
  pressure_(dd*pow(c0_,2),
	    wavelength,0,p0),
  xvelocity_(c0_*dd/d0,
	     wavelength,0,0),
  yvelocity_(0) {}

SpatialDistribution1D const& AcousticInitCond::getProfile
(string const& pname) const
{
  if("density"==pname)
    return density_;
  else if("pressure"==pname)
    return pressure_;
  else if("xvelocity"==pname)
    return xvelocity_;
  else if("yvelocity"==pname)
    return yvelocity_;
  else
    throw UniversalError("Unknown profile name "+pname);
}
