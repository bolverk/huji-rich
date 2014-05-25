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

  //! \brief Density
  double Density;

  //! \brief Pressure
  double Pressure;

  //! \brief Thermal energy per unit mass
  double Energy;

  //! \brief Velocity
  Vector2D Velocity;

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

Primitive operator+(Primitive const& p1, Primitive const& p2);

Primitive operator-(Primitive const& p1, Primitive const& p2);

Primitive operator/(Primitive const& p, double s);

Primitive operator*(Primitive const& p, double s);

Primitive operator*(double s, Primitive const& p);

Conserved operator*(double d, Conserved const& c);

Conserved operator+(Conserved const& c1, Conserved const& c2);

Conserved operator-(Conserved const& c1, Conserved const& c2);

Conserved operator/(Conserved const& c, double d);

double TotalEnergyDensity(Primitive const& p);

Conserved Primitive2Conserved(Primitive const& p);

Conserved Primitive2Conserved(Primitive const& p, double volume);

Conserved Primitive2Flux(Primitive const& p, Vector2D const& n);

#endif // HYDRODYNAMIC_VARIABLES_HPP
