#include <cmath>
#include "simple_waves_ideal_gas.hpp"

double calc_entropy(double d, double p, double g)
{
  return p/pow(d,g);
}

ConstEntropy::ConstEntropy(SpatialDistribution1D const& density,
			   double s, double g):
  density_(density), s_(s), g_(g) {}

double ConstEntropy::operator()(double x) const
{
  const double d = density_(x);
  return s_*pow(d,g_);
}

SoundSpeedDist::SoundSpeedDist
(SpatialDistribution1D const& pressure,
 SpatialDistribution1D const& density,
 EquationOfState const& eos):
  pressure_(pressure), density_(density), eos_(eos) {}

double SoundSpeedDist::operator()(double x) const
{
  const double p = pressure_(x);
  const double d = density_(x);
  return eos_.dp2c(d,p);
}

double calc_riemann_invariant(double v, double c,
			      double g, bool dir)
{
  double fdir = 2*dir-1;
  return v+fdir*2*c/(g-1);
}

ConstRiemannInv::ConstRiemannInv(double rv, bool dir,
				 SpatialDistribution1D const& sound_speed,
				 double g):
  rv_(rv), dir_(dir), sound_speed_(sound_speed), g_(g) {}

double ConstRiemannInv::operator()(double x) const
{
  const double c = sound_speed_(x);
  return rv_-calc_riemann_invariant(0,c,g_,dir_);
}

SimpleWaveIdealGasInitCond::SimpleWaveIdealGasInitCond
(SpatialDistribution1D const& density,
 double entropy, double adiabatic_index,
 double edge):
  density_(density),
  pressure_(density_,entropy,adiabatic_index),
  eos_(adiabatic_index),
  sound_speed_(pressure_,
	       density_,
	       eos_),
  c_edge_(sound_speed_(edge)),
  rv_edge_(calc_riemann_invariant(0,c_edge_,adiabatic_index,0)),
  xvelocity_(rv_edge_,0,sound_speed_,adiabatic_index),
  yvelocity_(0) {}

IdealGas const& SimpleWaveIdealGasInitCond::getEOS(void) const
{
  return eos_;
}

SpatialDistribution1D const& SimpleWaveIdealGasInitCond::getProfile(string pname) const
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
    throw UniversalError("Unknown variable name "+pname);
}

EntropyProf::EntropyProf(SpatialDistribution1D const& density,
			 SpatialDistribution1D const& pressure,
			 double adiabatic_index):
  density_(density), pressure_(pressure),
  g_(adiabatic_index) {}

double EntropyProf::operator()(double x) const
{
  const double d = density_(x);
  const double p = pressure_(x);
  return calc_entropy(d,p,g_);
}
