#include "sine_wave.hpp"
#include <cmath>

SineWave::SineWave(double amplitude,
		   double wavelength,
		   double phase,
		   double offset):
  amp_(amplitude),
  k_(2*M_PI/wavelength),
  ph_(phase),
  offset_(offset) {}

double SineWave::operator()(double x) const
{
  return amp_*sin(k_*x+ph_)+offset_;
}
