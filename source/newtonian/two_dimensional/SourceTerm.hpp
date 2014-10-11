/*! \file SourceTerm.hpp
  \brief Abstract class for source terms
  \author Elad Steinberg
*/

#ifndef SOURCETERM_HPP
#define SOURCETERM_HPP 1
#include "../../tessellation/tessellation.hpp"
#include "../common/hydrodynamic_variables.hpp"
#include "HydroBoundaryConditions.hpp"
#include "../../misc/utils.hpp"
#include "physical_geometry.hpp"

//! \brief Abstract class for external forces
class SourceTerm
{
public:
  /*!
    \brief Calcualtes the change in conserved variables done on a cell from a source term
    \param tess The tessellation
    \param pg Physical geometry
    \param cells The hydrodynmic variables of the cell
    \param point The index of the cell's CM
    \param fluxes Hydrodynmic fluxes
    \param point_velocity Velocities of mesh generating points
    \param hbc Hydrodynamic boundary condition
    \param tracers The intensive scalar tracers
	\param dtracer The change per unit time of the extensive tracers, doesn't matter the input since the result overwrites the vector
	\param lengthes THe corrected length of the edges
	\param t Time
    \param dt The time step size
    \return The flux of conserved variables
  */
  virtual Conserved Calculate
  (Tessellation const& tess,
   const PhysicalGeometry& pg,
   vector<Primitive> const& cells,
   int point,
   vector<Conserved> const& fluxes,
   vector<Vector2D> const& point_velocity,
   HydroBoundaryConditions const& hbc,
   vector<vector<double> > const& tracers,
   vector<double>& dtracer,vector<double> const& lengthes,
   double t,
   double dt)=0;

  virtual ~SourceTerm(void);
};

#endif //SOURCETERM_HPP
