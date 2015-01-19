#ifndef CONSERVED_3D_HPP
#define CONSERVED_3D_HPP 1

#include "../../3D/GeometryCommon/Vector3D.hpp"

//! \brief Conserved variables for a 3D computational cell
class Conserved3D
{
public:

  //! \brief Mass
  double mass;

  //! \brief Momentum
  Vector3D momentum;

  //! \brief Energy
  double energy;

  //! \brief Tracers
  vector<double> tracers;

  //! \brief Class constructor (sets everything to zero)
  Conserved3D(void);

  /*! \brief Class constructor (does not initialize tracers)
    \param mass_i Mass
    \param momentum_i Momentum
    \param energy_i Energy
   */
  Conserved3D(double mass_i,
	      const Vector3D& momentum_i,
	      double energy_i);

  /*! \brief Class constructor
    \param mass_i Mass
    \param momentum_i Momentum
    \param energy_i Energy
    \param tracers_i Tracers
   */
  Conserved3D(double mass_i,
	      const Vector3D& momentum_i,
	      double energy_i,
	      const vector<double>& tracers_i);

  /*! \brief Reduction operator
    \param diff Difference
    \return Reference to self
   */
  Conserved3D& operator-=(const Conserved3D& diff);

  /*! \brief Addition operator
    \param diff Difference
    \return Reference to self
   */
  Conserved3D& operator+=(const Conserved3D& diff);
};

/*! \brief Scalar product operator
  \param s Scalar
  \param c Conserved variable
  \return Product of scalar with conserved
 */
Conserved3D operator*(double s, const Conserved3D& c);

/*! \brief Scalar division operator
  \param c Conserved variable
  \param s Scalar
  \return Ratio
 */
Conserved3D operator/(const Conserved3D& c, double s);

#endif // CONSERVED_3D_HPP
