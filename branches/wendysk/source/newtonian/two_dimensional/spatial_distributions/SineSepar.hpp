/*! \brief A spatial distribution with a sine wave function as a separatrix
  \author Elad Steinberg
 */

#ifndef SINESEPAR_HPP
#define SINESEPAR_HPP 1

#include "../spatial_distribution2d.hpp"
#include <cmath>

class SineSepar: public SpatialDistribution
{
private:

  double Amplitude,freq,offset,above_,under_;

public:

  /*! \brief Class constructor
    \param Amp Amplitude
    \param f Frequency
	\param off Offset
    \param above Value above the sine
    \param under Value under the sine
   */
  SineSepar(double Amp,double f,double off,double above,double under);
  ~SineSepar();
  double EvalAt(Vector2D const& r) const;
};

#endif // SINESEPAR_HPP
