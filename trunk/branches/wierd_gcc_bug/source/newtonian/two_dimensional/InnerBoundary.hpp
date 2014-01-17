//Pure virtual
#ifndef INNERBOUNDARY_CONDITIONS_HPP
#define INNERBOUNDARY_CONDITIONS_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../common/riemann_solver.hpp"
#include "../../tessellation/tessellation.hpp"
#include "HydroBoundaryConditions.hpp"
#include <cmath>

/*! \brief Inner Hydro Boundary Conditions
  \author Elad Steinberg
 */
class InnerBoundaryConditions: public HydroBoundaryConditions
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
	\return The flux
   */
  virtual Conserved CalcFlux
  (Tessellation const* tessellation,
   vector<Primitive> const& cells,Vector2D const& edge_velocity,
   Edge const& edge,
   SpatialReconstruction const* interp,double dt,double time) const=0;

  /*! \brief Calculates the velocity of the boundary edge
    \param tessellation Point and edge positions
	\param point_velocities Velocities of the mesh generating points
    \param edge Boundary edge    
	\param time The sim time
	\return Velocity of edge
   */
  virtual Vector2D CalcEdgeVelocity(Tessellation const* tessellation,
				    vector<Vector2D> const& point_velocities,
				    Edge const& edge,double time) const = 0;
  /*!
  \brief Checks if the edge is on a boundary
  \param Data Point and edge positions
   \param edge The edge to check
   \returns If it is on a boundary or not
   */
  virtual bool IsBoundary(Edge const& edge,Tessellation const* Data)const=0;
  /*!
  \brief Checks if the cell is a ghost cell
  \param i The index of the cell
  \param Data Point and edge positions
   \returns If the cell is a ghost cell or not
   */
  virtual bool IsGhostCell(int i,Tessellation const* Data) const=0;

  /*!
  \brief Returns the number of inner points
  \return The number of inner points
  */
  virtual int GetPointNum(void)const=0;

  virtual Primitive GetBoundaryPrimitive(Edge const& edge,
	  Tessellation const* Data,vector<Primitive> const& cells,
	  double time)const=0;
  /*!
  \brief Returns the tracers on the other side of the boundary edge
  \param edge The edge to check
  \param Data The tessellation
  \param tracers The tracers
  \param time The sim time
  \return The tracers on the other side
  */
  virtual vector<double> GetBoundaryTracers(Edge const& edge,
	  Tessellation const* Data,vector<vector<double> > const& tracers,
	  double time)const=0;

  /*! \brief Calculates the flux on the boundary edge
    \param tessellation Point and edge positions
    \param tracers The tracers
	\param dm The mass flux through the edge
    \param edge Boundary edge
	\param index Tracer index
	\param dt The time step
	\param time The sim time
	\param interp Scalar interpolation
	\return The flux of the tracer
   */

  vector<double> CalcTracerFlux
  (Tessellation const* tessellation,
   vector<vector<double> > const& tracers,double dm,
   Edge const& edge,int index,double dt,
   double time,
   SpatialReconstruction const* interp) const;
};

#endif // INNERBOUNDARY_CONDITIONS_HPP
