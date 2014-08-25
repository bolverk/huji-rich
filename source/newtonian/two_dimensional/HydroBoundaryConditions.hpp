/*! \file HydroBoundaryConditions.hpp
  \brief Hydro Boundary Conditions
  \author Elad Steinberg
 */

#ifndef HYDRO_BOUNDARY_CONDITIONS_HPP
#define HYDRO_BOUNDARY_CONDITIONS_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../common/riemann_solver.hpp"
#include "../../tessellation/tessellation.hpp"
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
    \param edge_velocity Velocity of the edge
    \param edge Boundary edge
	\param interp The spatial reconstruction
	\param dt The time step
	\param time The sim time
	\return The flux through the edge
   */
  virtual Conserved CalcFlux(Tessellation const& tessellation,
	  vector<Primitive> const& cells,Vector2D const& edge_velocity,
	  Edge const& edge,SpatialReconstruction const& interp,double dt,
	  double time) const = 0;

  /*! \brief Calculates the velocity of the boundary edge
    \param tessellation Point and edge positions
	\param point_velocities Velocities of the mesh generating points
    \param edge Boundary edge
    \param time The sim time
	\return Velocity of edge
   */
  virtual Vector2D CalcEdgeVelocity
  (Tessellation const& tessellation,
   vector<Vector2D> const& point_velocities,
   Edge const& edge,
   double time) const = 0;

  /*!
  \brief Checks if the edge is on a boundary
  \param tess The tessellation
   \param edge The edge to check
   \returns If it is on a boundary or not
   */
  virtual bool IsBoundary(Edge const& edge,Tessellation const& tess)const=0;

  /*!
  \brief Checks if the cell is a ghost cell
  \param i The index of the cell
  \param tess The tessellation
   \returns If the cell is a ghost cell or not
   */
  virtual bool IsGhostCell(int i,Tessellation const& tess) const=0;
  /*!
  \brief Returns the normal direction to an edge
  \param edge The edge
  \param tess The tessellation
   \returns The normal direction
   */
  Vector2D Normal(Edge const& edge, Tessellation const& tess)const;

  /*!
  \brief Returns the primitive on the other side of the boundary edge
  \param edge The edge to check
  \param tess The tessellation
  \param cells The primitive cells
  \param time The sim time
  \return The primitive on the other side
  */
  virtual Primitive GetBoundaryPrimitive(Edge const& edge,
	  Tessellation const& tess,vector<Primitive> const& cells,
	  double time)const=0;
  /*!
  \brief Returns the tracers on the other side of the boundary edge
  \param edge The edge to check
  \param tess The tessellation
  \param tracers The tracers
  \param time The sim time
  \return The tracers on the other side
  */
  virtual vector<double> GetBoundaryTracers(Edge const& edge,
	  Tessellation const& tess,vector<vector<double> > const& tracers,
	  double time)const=0;

  /*! \brief Calculates the tracer flux on the boundary edge
    \param tessellation Point and edge positions
	\param cells The hydro cells
    \param tracers The tracers
	\param dm The mass flux through the edge
    \param edge Boundary edge
	\param index The cell's index
	\param dt The time step
	\param time The sim time
	\param interp Scalar interpolation method
	\param edge_velocity The velocity of the edge
	\return The flux of the tracer
   */
  virtual vector<double> CalcTracerFlux
  (Tessellation const& tessellation,
   vector<Primitive> const& cells,
   vector<vector<double> > const& tracers,
   double dm,
   Edge const& edge,
   int index,
   double dt,
   double time,
   SpatialReconstruction const& interp,
   Vector2D const& edge_velocity
   ) const = 0;

  //! \brief virtual destructor
  virtual ~HydroBoundaryConditions(void);
};

#endif // HYDRO_BOUNDARY_CONDITIONS_HPP