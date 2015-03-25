/*! \file simple_waves_ideal_gas.hpp
  \brief Initial conditions for simple waves problem in an ideal gas
  \author Almog Yalinewich
 */

#ifndef SIMPLE_WAVES_IDEAL_GAS_HPP
#define SIMPLE_WAVES_IDEAL_GAS_HPP 1

#include <string>
#include "../one_dimensional/spatial_distribution1d.hpp"
#include "../common/equation_of_state.hpp"
#include "../common/ideal_gas.hpp"

//using namespace std;
using std::string;

/*! \brief Calculates the entropy of an ideal gas
  \param d Density
  \param p Pressure
  \param g Adiabatic index
  \return Entropy
 */
double calc_entropy(double d, double p, double g);

/*! \brief Pressure distribution with constant entropy construct
 */
class ConstEntropy: public SpatialDistribution1D
{
public:

  /*! \brief Class constructor
    \param density Density distribution
    \param s Entropy
    \param g Adiabatic index
   */
  ConstEntropy(SpatialDistribution1D const& density,
	       double s, double g);

  double operator()(double x) const;

private:

  SpatialDistribution1D const& density_;
  const double s_;
  const double g_;
};

/*! \brief Spatial distribution of the speed of sound
 */
class SoundSpeedDist: public SpatialDistribution1D
{
public:

  /*! \brief Class constructor
    \param pressure Pressure Distribution
    \param density Density distribution
    \param eos Equation of state
   */
  SoundSpeedDist(SpatialDistribution1D const& pressure,
		 SpatialDistribution1D const& density,
		 EquationOfState const& eos);

  double operator()(double x) const;

private:

  SpatialDistribution1D const& pressure_;
  SpatialDistribution1D const& density_;
  EquationOfState const& eos_;
};

/*! \brief Calculates the Riemann invariant
  \param v Velocity
  \param c Sound speed
  \param g Adiabatic index
  \param dir Direction
  \return Riemann invariant
 */
double calc_riemann_invariant(double v, double c,
			      double g, bool dir);

//! \brief Velocity distribution where the Riemann invariant is constant throughout
class ConstRiemannInv: public SpatialDistribution1D
{
public:

  /*! \brief Class constructor
    \param rv Value of the riemann invariant
    \param dir Direction
    \param sound_speed Sound speed distribution
    \param g Adiabatic index
   */
  ConstRiemannInv(double rv, bool dir,
		  SpatialDistribution1D const& sound_speed,
		  double g);

  double operator()(double x) const;

private:

  const double rv_;
  const bool dir_;
  SpatialDistribution1D const& sound_speed_;
  const double g_;
};

/*! \brief Initial conditions for simple waves in an ideal gas
  \details The user supplies this class with a density distribution, and the class calculates the corresponding pressure and velocity so it would be a simple wave.
 */
class SimpleWaveIdealGasInitCond
{
public:

  /*! \brief Class constructor
    \param density Density profile
    \param entropy Entropy
    \param adiabatic_index Adiabatic index
    \param edge Position outside the disturbed region
   */
  SimpleWaveIdealGasInitCond
  (SpatialDistribution1D const& density,
   double entropy,
   double adiabatic_index,
   double edge=3);

  /*! \brief Returns the equation of state
   */
  IdealGas const& getEOS(void) const;

  /*! \brief Returns the appropriate spatial profile
    \param pname Name of the Hydrodynamic variables
    \return Spatial profile
   */
  SpatialDistribution1D const& getProfile(string pname) const;

private:

  SpatialDistribution1D const& density_;
  ConstEntropy pressure_;
  IdealGas eos_;
  SoundSpeedDist sound_speed_;
  double c_edge_;
  double rv_edge_;
  ConstRiemannInv xvelocity_;
  Uniform yvelocity_;
};

//! \brief Spatial distribution of the entropy
class EntropyProf: public SpatialDistribution1D
{
public:

  /*! \brief Class constructor
    \param density Density distribution
    \param pressure Pressure distribution
    \param adiabatic_index Adiabatic index
   */
  EntropyProf(SpatialDistribution1D const& density,
	      SpatialDistribution1D const& pressure,
	      double adiabatic_index);

  double operator()(double x) const;

private:

  SpatialDistribution1D const& density_;
  SpatialDistribution1D const& pressure_;
  double g_;
};

#endif // SIMPLE_WAVES_IDEAL_GAS_HPP
