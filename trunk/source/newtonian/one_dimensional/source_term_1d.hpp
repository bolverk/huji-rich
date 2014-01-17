/*! \file source_term_1d.hpp
  \brief Abstract class for external forces
  \author Almog Yalinewich
*/

#ifndef EXTERNAL_FORCES_1D_HPP
#define EXTERNAL_FORCES_1D_HPP 1

#include <vector>
#include "../common/hydrodynamic_variables.hpp"

using namespace std;

//! \brief Abstract class for external forces
class ExternalForces1D
{
public:

  /*! \brief Calculates the change in the extensive conserved variables due to external forces
    \param vertices Position of the vertices
    \param cells Hydrodynamic cells
    \param point Point in which the forces will be evaluated
    \param t Simulation time
    \param dt Time step
    \return Value of the external force
   */
  virtual Conserved calc
  (vector<double> const& vertices,
   vector<Primitive> const& cells,
   int point,
   double t,
   double dt) const = 0;

  virtual ~ExternalForces1D(void);
};

#endif // EXTERNAL_FORCES_1D_HPP
