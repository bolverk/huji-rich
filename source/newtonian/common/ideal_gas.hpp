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

  double g_;

public:

  /*! \brief Class constructor
    \param AdiabaticIndex Adiabatic index
   */
  explicit IdealGas(double AdiabaticIndex);

  /*! \brief Returns the adiabatic index
    \return Adiabatic index
   */
  double getAdiabaticIndex(void) const;

  double dp2e(double d, double p, const boost::container::flat_map<string,double>& tracers) const;

  double de2p(double d, double e, const boost::container::flat_map<string,double>& tracers) const;

  double dp2c(double d, double p, const boost::container::flat_map<string,double>& tracers) const;

  double de2c(double d, double e, const boost::container::flat_map<string,double>& tracers) const;

  double dp2s(double d, double p, const boost::container::flat_map<string,double>& tracers) const;

  double sd2p(double s, double d, const boost::container::flat_map<string,double>& tracers) const;
};

#endif // IDEAL_GAS_HPP
