//! \file hydrodynamic_variables.hpp
//! \brief Hydrodynamic variables
//! \author Almog Yalinewich

#ifndef HYDRODYNAMIC_VARIABLES_HPP
#define HYDRODYNAMIC_VARIABLES_HPP 1

#include "../../tessellation/geometry.hpp"

//! \brief Set of conserved variables (extensive)
class Conserved
{
public:

  //! \brief Null constructor (sets all members to zero)
  Conserved(void);

  /*! \brief Class constructor
    \param mass Mass
    \param momentum Momentum
    \param energy Energy
   */
  Conserved(double mass,
	    Vector2D const& momentum,
	    double energy);

  Conserved(const Conserved& other);

  //! \brief Mass
  double Mass;

  //! \brief Momentum
  Vector2D Momentum;

  //! \brief Total energy (kinetic + thermal)
  double Energy;

  /*! \brief Addition of flux
    \param f Flux
    \return Sum
   */
  Conserved& operator+=(Conserved const& f);

  /*! \brief subtraction of flux
    \param c Conserved variables
    \return Difference
   */
  Conserved& operator-=(Conserved const& c);

  /*! \brief Assigmnet operator
	\param other The object to be copied
	\return The copy
  */
  Conserved& operator=(Conserved const&other);
};

//! \brief Primitive hydrodynamic variables
class Primitive
{
public:

  enum {NUMBER_PRIMITIVE_VARIABLES=6};

  Primitive(void);

  /*! \brief Class constructor
    \param density_i Density
    \param pressure_i Pressure
    \param velocity_i Velocity
    \param energy_i Energy
    \param sound_speed_i Sound speed
   */
  Primitive(double density_i,
	    double pressure_i,
	    Vector2D const& velocity_i,
	    double energy_i,
	    double sound_speed_i);

  Primitive(const Primitive& other);

  //! \brief Density
  double Density;

  //! \brief Pressure
  double Pressure;

  //! \brief Velocity
  Vector2D Velocity;

  //! \brief Thermal energy per unit mass
  double Energy;

  //! \brief Speed of sound
  double SoundSpeed;

  //! \brief Returns the numbers of members
  int GetVarNo(void) const;

  /*! \brief Operator for adding members of two primitives
    \param p Another primitive
    \return Reference to primitive variables
   */
  Primitive& operator+=(Primitive const& p);

  /*! \brief Subscript operator
    \details Enumerates the variables
    0 - Density
    1 - Pressure
    2 - Energy
    3 - Sound speed
    4 - x Velocity
    5 - y Velocity
    \param index Parameter index
    \return Reference to variable
   */
  double& operator[](int index);

  /*! \brief Subscript operator
    \details Enumerates the variables
    0 - Density
    1 - Pressure
    2 - Energy
    3 - Sound speed
    4 - x Velocity
    5 - y Velocity
    \param index Parameter index
    \return Value of variable
   */
  double operator[](int index) const;

  /*! \brief Assigmnet operator
	\param other The object to be copied
	\return The copy
  */
  Primitive& operator=(Primitive const&other);
};

/*! \brief Checks if on of the fields of Primitive is a nan
  \param p Primitive variables
  \return True if a primitive has a nan
 */
bool primitive_has_nan(Primitive const& p);

/*! \brief Term by term addition
  \param p1 Primitive variables
  \param p2 Primitive variables
  \return Primitive variables
 */
Primitive operator+(Primitive const& p1, Primitive const& p2);

/*! \brief Term by term subtraction
  \param p1 Primitive variables
  \param p2 Primitive variables
  \return Primitive variables
 */
Primitive operator-(Primitive const& p1, Primitive const& p2);

/*! \brief Scalar division
  \param p Primitive
  \param s Scalar
  \return Primitive variables
 */
Primitive operator/(Primitive const& p, double s);

/*! \brief Scalar multiplication on the right
  \param p Primitive
  \param s Scalar
  \return Primitive variables
 */
Primitive operator*(Primitive const& p, double s);

/*! \brief Scalar multiplication on the left
  \param s Scalar
  \param p Primitive variables
  \return Primitive variables
 */
Primitive operator*(double s, Primitive const& p);

/*! \brief Scalar multiplication
  \param d Scalar
  \param c Conserved variables
  \return Conserved variables
 */
Conserved operator*(double d, Conserved const& c);

/*! \brief Term by term addition
  \param c1 Conserved variables
  \param c2 Conserved variables
  \return Conserved variables
 */
Conserved operator+(Conserved const& c1, Conserved const& c2);

/*! \brief Term by term subtraction
  \param c1 Conserved variables
  \param c2 Conserved variables
  \return Conserved variables
 */
Conserved operator-(Conserved const& c1, Conserved const& c2);

/*! \brief Division by a scalar
  \param c Conserved variables
  \param d Scalar
  \return Conserved variables
 */
Conserved operator/(Conserved const& c, double d);

/*! \brief Calculates the total energy density
  \param p Primitive variables
  \return Total energy density
 */
double TotalEnergyDensity(Primitive const& p);

/*! \brief Converts primitive variables to conserved intensive
  \param p Primitive variables
  \return Conserved variables
 */
Conserved Primitive2Conserved(Primitive const& p);

/*! \brief Converts primitive variables to conserved extensive
  \param p Primitive variables
  \param volume Cell volume
  \return Conserved variables
 */
Conserved Primitive2Conserved(Primitive const& p, double volume);

/*! \brief Converts primitive variables to flux
  \param p Primitive variables
  \param n Normal to surface through which the flux passes
  \return Conserved variables
 */
Conserved Primitive2Flux(Primitive const& p, Vector2D const& n);

#endif // HYDRODYNAMIC_VARIABLES_HPP
