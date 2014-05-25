#ifndef RIGIDWALLHYDRO_HPP
#define RIGIDWALLHYDRO_HPP 1

#include "../HydroBoundaryConditions.hpp"
#include "../InnerBoundary.hpp"
#include "../hydrodynamics_2d.hpp"

/*! \brief Rigid Wall Hydro Boundary Conditions
\author Elad Steinberg
*/
class RigidWallHydro: public HydroBoundaryConditions
{
public:
	/*! \brief Class constructor
	\param rs The Riemann solver
	*/
  RigidWallHydro(RiemannSolver const& rs);
  //! \brief Class destructor
  ~RigidWallHydro();

  Conserved CalcFlux(Tessellation const* tessellation,
	  vector<Primitive> const& cells,Vector2D const& edge_velocity,
	  Edge const& edge,SpatialReconstruction const* interp,double dt,
	  double time) const;

  Vector2D CalcEdgeVelocity
  (Tessellation const* tessellation,
   vector<Vector2D> const& point_velocities,
   Edge const& edge,
   double time) const;

  bool IsBoundary(Edge const& edge,Tessellation const* Data)const;

  bool IsGhostCell(int i,Tessellation const* Data) const;
  /*!
  \brief Returns the real cell on a boundary edge
  \param edge The edge
  \return The index of the real cell
  */
  int GetRealCell(Edge const& edge)const;
  /*!
	\brief Calcualtes the flux between the real cell and the rigid counterpart
	\param tessellation The tessellation
	\param cells The primitive cells
	\param edge_velocity The velocities of the edges
	\param rs The Riemann solver
	\param edge The edge to calculate
	\param interp The spatial interpolation
	\param dt The time step
	\param ci The index of the real cell in teh edge
	\return The flux
	*/
  Conserved CalcFluxCi
  (Tessellation const* tessellation,
   vector<Primitive> const& cells,Vector2D const& edge_velocity,
   RiemannSolver const* rs,Edge const& edge,
   SpatialReconstruction const* interp,double dt,int ci) const;

  Primitive GetBoundaryPrimitive
  (Edge const& edge,
   Tessellation const* Data,
   vector<Primitive> const& cells,double time)const;

  vector<double> GetBoundaryTracers(Edge const& edge,
   Tessellation const* Data,
   vector<vector<double> > const& tracers,double time)const;
  
  vector<double> CalcTracerFlux
  (Tessellation const* tessellation,vector<Primitive> const& cells,
   vector<vector<double> > const& tracers,double dm,
   Edge const& edge,int index,double dt,
   double time,SpatialReconstruction const* interp) const;

private:
  
  RiemannSolver const& rs_;
};

#endif // RIGIDWALLHYDRO_HPP
