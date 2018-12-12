/*! \brief Base class for equation of state
  \file equation_of_state.hpp
  \author Almog Yalinewich
 */

#ifndef EQUATION_OF_STATE_HPP
#define EQUATION_OF_STATE_HPP 1

#include <string>
#include <array>
#include "../two_dimensional/computational_cell_2d.hpp"

/** typedef for tracer vector */
//typedef std::array<double,MAX_TRACERS> tvector;
/** typedef for string vector */
//typedef std::array<bool, MAX_STICKERS> svector;
using std::string;
using std::vector;

//! \brief Base class for equation of state
class EquationOfState
{
public:

  /*! \brief Calculates the thermal energy per unit mass
    \param d Density
    \param p Pressure
    \param tracers Tracers
	\param tracernames The names of the tracers
    \return Thermal energy per unit mass
   */
  virtual double dp2e(double d, double p, tvector const& tracers = tvector(),vector<string> const& tracernames = vector<string>())
	  const = 0;

  /*! \brief Calculates the pressure
    \param d Density
    \param e Specific thermal energy
    \param tracers Tracers
	\param tracernames The names of the tracers
    \return Presusre
   */
  virtual double de2p(double d, double e,
	  tvector const& tracers = tvector(), vector<string> const& tracernames = vector<string>()) const = 0;

  /*! \brief Calculates the speed of sound
    \param d Density
    \param e Specific thermal energy
    \param tracers Tracers
	\param tracernames The names of the tracers
    \return Speed of sound
  */
  virtual double de2c(double d, double e,
	  tvector const& tracers = tvector(), vector<string> const& tracernames = vector<string>()) const = 0;

  /*! \brief Calculates the speed of sound
    \param d Density
    \param p Pressure
    \param tracers Tracers
	\param tracernames The names of the tracers
    \return Speed of sound
   */
  virtual double dp2c(double d, double p,
	  tvector const& tracers = tvector(), vector<string> const& tracernames = vector<string>()) const = 0;

  /*! \brief Calculates the entropy per unit mass
    \param d Density
    \param p Pressure
    \param tracers Tracers
	\param tracernames The names of the tracers
    \return Entropy
  */
  virtual double dp2s(double d, double p,
	  tvector const& tracers = tvector(), vector<string> const& tracernames = vector<string>()) const = 0;

  /*! \brief Calculates the pressure from the netropy
    \param s Entropy
    \param d Density
    \param tracers Tracers
	\param tracernames The names of the tracers
    \return Entropy
  */
  virtual double sd2p(double s, double d,
	  tvector const& tracers = tvector(), vector<string> const& tracernames = vector<string>()) const = 0;

  virtual ~EquationOfState(void);
};

#endif
