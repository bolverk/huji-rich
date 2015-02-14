/*! \file computational_cell.hpp
  \brief A container for the hydrodynamic variables
  \author Almog Yalinewich
 */

#ifndef COMPUTATIONAL_CELL_HPP
#define COMPUTATIONAL_CELL_HPP 1

#include "../../3D/GeometryCommon/Vector3D.hpp"

//! \brief Container for the hydrodynamic variables
class ComputationalCell
{
public:

  //! \brief Density
  double density;

  //! \brief Pressure
  double pressure;

  //! \brief Velocity
  Vector3D velocity;

  //! \brief Tracers
  vector<double> tracers;

  //! \brief Null constructor
  ComputationalCell(void);

  /*! \brief Class constructor
    \param density_i Density
    \param pressure_i Pressure
    \param velocity_i Velocity
   */
  ComputationalCell(double density_i,
		    double pressure_i,
		    const Vector3D& velocity_i);

  /*! \brief Class constructor
    \param density_i Density
    \param pressure_i Pressure
    \param velocity_i Velocity
    \param tracers_i Tracers
   */
  ComputationalCell(double density_i,
		    double pressure_i,
		    const Vector3D& velocity_i,
		    const vector<double>& tracers_i);
};

#endif // COMPUTATIONAL_CELL_HPP
