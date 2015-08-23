/*! \file HydroBoundaryConditions.hpp
  \brief Hydro Boundary Conditions
  \author Elad Steinberg
*/

#ifndef HYDRO_BOUNDARY_CONDITIONS_HPP
#define HYDRO_BOUNDARY_CONDITIONS_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../common/riemann_solver.hpp"
#include "../../tessellation/tessellation.hpp"
#include "extensive.hpp"
#include "spatial_reconstruction.hpp"
#include <cmath>
#include <algorithm>
#include <functional>

//! \brief Square box outer boundary conditions with two sides reflective and two periodic. The x direction is taken to be periodic.
class HydroBoundaryConditions
{
public:

  /*! \brief Calculates the flux on the boundary edge
    \param tessellation Point and edge positions
    \param cells Hydrodynamic variables
    \return List of indices and corresponding fluxes
  */
  virtual vector<pair<size_t,Extensive> > operator()
  (const Tessellation& tessellation,
   const vector<ComputationalCell>& cells) const = 0;

  //! \brief virtual destructor
  virtual ~HydroBoundaryConditions(void);
};

#endif // HYDRO_BOUNDARY_CONDITIONS_HPP
