/*! \file collela.hpp
  \brief Smooth spatial distribution
  \author Almog Yalinewich
 */

#ifndef COLLELA_HPP
#define COLLELA_HPP 1

#include "../one_dimensional/spatial_distribution1d.hpp"

/*! \brief A smooth spatial distribution attributed to Philip Collela
  \details For more information refer to section 4.6 in <a href="http://adsabs.harvard.edu/abs/2006ApJS..164..255Z"> W. Zhang, A. MacFadyen, "RAM: A Relativistic Adaptive Mesh Refinement Hydrodynamics Code", ApJS 164, 255 (2006) </a>.
 */
class Collela: public SpatialDistribution1D
{
public:
	/*!
	\brief Class constructor
	\param ref Unperturbed value
	\param a Relative amplitude of the wave
	\param l Wavelength
	\param offset Horizontal offset
	*/
  Collela(double ref, double a, double l, double offset);

  double operator()(double x) const;

private:

  const double ref_;
  const double a_;
  const double l_;
  const double offset_;
};

#endif // COLLELA_HPP
