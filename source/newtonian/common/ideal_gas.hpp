/*! \file ideal_gas.hpp
  \brief Ideal gas equation of state
  \author Almog Yalinewich
 */

#ifndef IDEAL_GAS_HPP
#define IDEAL_GAS_HPP 1

#include "equation_of_state.hpp"
#include "../../misc/utils.hpp"

//! \brief Ideal gas equation of state
class IdealGas: public EquationOfState
{
private:

  double const g_, f_, beta_, mu_;

public:

  /*! \brief Class constructor
    \param AdiabaticIndex Adiabatic index
   */
  explicit IdealGas(double const AdiabaticIndex);

  /*! \brief Class constructor for ideal gas with temperature equation as well, e=f*T^beta*rho^{-mu}
    \param AdiabaticIndex Adiabatic index
    \param f The constant prefactor
    \param beta The power law of the temperature
    \param mu The inverse of the power law of the density
   */
  IdealGas(double const AdiabaticIndex, double const f, double const beta, double const mu);

  /*! \brief Returns the adiabatic index
    \return Adiabatic index
   */
  double getAdiabaticIndex(void) const;

  double dp2e(double d, double p, tvector const& tracers, vector<string> const& tracernames) const override;

  double de2p(double d, double e, tvector const& tracers, vector<string> const& tracernames) const override;

  double dp2c(double d, double p, tvector const& tracers, vector<string> const& tracernames) const override;

  double de2c(double d, double e, tvector const& tracers, vector<string> const& tracernames) const override;

  double dp2s(double d, double p, tvector const& tracers, vector<string> const& tracernames) const override;

  double sd2p(double s, double d, tvector const& tracers, vector<string> const& tracernames) const override;

  double dT2cv(double const d, double const T, tvector const& tracers, vector<string> const& tracernames) const override;

  double de2T(double const d, double const T, tvector const& tracers, vector<string> const& tracernames) const override;

  double dT2e(double const d, double const T, tvector const& tracers, vector<string> const& tracernames) const override;
};

#endif // IDEAL_GAS_HPP
