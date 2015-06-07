/*! \file sine_wave.hpp
  \brief Single mode sine wave spatial distribution
  \author Almog Yalinewich
 */

#ifndef SINE_WAVE_HPP
#define SINE_WAVE_HPP 1

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif // _MSC_VER
#include "../one_dimensional/spatial_distribution1d.hpp"
//! \brief Sine wave spatial distribution
class SineWave: public SpatialDistribution1D
{
public:

  /*!
    \brief Class constructor return amp*sin(x*k+phase)+offest
    \param amplitude The amplitude of the sine function
    \param wavelength The wavelength
    \param phase The phase
    \param offset The constant added to the sine
  */
  SineWave(double amplitude,
	   double wavelength,
	   double phase,
	   double offset);

  double operator()(double x) const;

private:

  double amp_;
  double k_;
  double ph_;
  double offset_;
};

#endif // SINE_WAVE_HPP
