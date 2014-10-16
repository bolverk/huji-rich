/*! \file acoustic.hpp
  \brief Initial conditions for single mode acoustic waves
  \author Almog Yalinewich
 */

#ifndef ACOUSTIC_HPP
#define ACOUSTIC_HPP 1

#include <string>
#include "../common/equation_of_state.hpp"
#include "../one_dimensional/spatial_distribution1d.hpp"
#include "sine_wave.hpp"

using std::string;

/*! \brief Initial conditions giving rise to acoustic waves
 */
class AcousticInitCond
{
public:

  /*! \brief Class constructor
    \param d0 Mean density
    \param p0 Mean pressure
    \param eos Equation of state
    \param dd Density perturbation amplitude
    \param wavelength Wavelength
   */
  AcousticInitCond(double d0, double p0,
		   EquationOfState const& eos,
		   double dd, double wavelength);

  /*! \brief Returns the spatial distribution of one of the hydrodynamic variables
    \param pname Name of the hydrodynamic variable
	\return The profile
   */
  SpatialDistribution1D const& getProfile(string const& pname) const;

private:

  double c0_;
  SineWave density_;
  SineWave pressure_;
  SineWave xvelocity_;
  Uniform yvelocity_;
};

#endif // ACOUSTIC_HPP
