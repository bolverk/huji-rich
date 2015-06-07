/*! \brief Base class for equation of state
  \file equation_of_state.hpp
  \author Almog Yalinewich
 */

#ifndef EQUATION_OF_STATE_HPP
#define EQUATION_OF_STATE_HPP 1

#include <map>
#include <string>

using std::map;
using std::string;

//! \brief Base class for equation of state
class EquationOfState
{
public:

  /*! \brief Calculates the thermal energy per unit mass
    \param d Density
    \param p Pressure
    \return Thermal energy per unit mass
   */
  virtual double dp2e(double d, double p, const map<string,double>& tracers=map<string,double>()) const = 0;

  /*! \brief Calculates the pressure
    \param d Density
    \param e Specific thermal energy
    \return Presusre
   */
  virtual double de2p(double d, double e,
		      const map<string,double>& tracers=map<string,double>()) const = 0;

  /*! \brief Calculates the speed of sound
    \param d Density
    \param e Specific thermal energy
    \return Speed of sound
  */
  virtual double de2c(double d, double e,
		      const map<string,double>& tracers=map<string,double>()) const = 0;

  /*! \brief Calculates the speed of sound
    \param d Density
    \param p Pressure
    \return Speed of sound
   */
  virtual double dp2c(double d, double p,
		      const map<string,double>& tracers=map<string,double>()) const = 0;

  /*! \brief Calculates the entropy per unit mass
    \param d Density
    \param p Pressure
    \return Entropy
  */
  virtual double dp2s(double d, double p,
		      const map<string,double>& tracers=map<string,double>()) const = 0;
  /*! \brief Calculates the pressure from the netropy
    \param s Entropy
    \param d Density
    \return Entropy
  */
  virtual double sd2p(double s, double d,
		      const map<string,double>& tracers=map<string,double>()) const = 0;

  virtual ~EquationOfState(void);
};

#endif
