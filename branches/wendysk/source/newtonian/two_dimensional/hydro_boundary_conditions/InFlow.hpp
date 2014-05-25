#ifndef INFLOW_HPP
#define INFLOW_HPP 1

#include "../HydroBoundaryConditions.hpp"
#include "../InnerBoundary.hpp"

/*! \brief Inflow Hydro Boundary Conditions
\author Elad Steinberg
*/
class InFlow: public HydroBoundaryConditions
{
public:
	/*! \brief Class constructor
	\param InFlux The primitive in to mimic on the outer side of the boundary
	\param rs The Riemann solver
	\param outer_tracer The tracer in to mimic on the outer side of the boundary
	*/
  InFlow(Primitive const& InFlux,RiemannSolver const& rs,
	vector<double> outer_tracer=vector<double>());
  //! \brief Class destructor
  ~InFlow();

  Conserved CalcFlux
  (Tessellation const* tessellation,
   vector<Primitive> const& cells,Vector2D const& edge_velocity,
   Edge const& edge,
   SpatialReconstruction const* interp,double dt,
   double time)const;

  Vector2D CalcEdgeVelocity(Tessellation const* tessellation,
			    vector<Vector2D> const& point_velocities,
			    Edge const& edge, double time) const;
  
  bool IsBoundary(Edge const& edge,Tessellation const* Data)const;

  bool IsGhostCell(int i,Tessellation const* Data) const;

  Primitive GetBoundaryPrimitive
  (Edge const& edge,
   Tessellation const* Data,
   vector<Primitive> const& cells,double time)const;

  vector<double> GetBoundaryTracers(Edge const& edge,Tessellation const* Data,
	vector<vector<double> > const& tracers,double time)const;

 vector<double> CalcTracerFlux(Tessellation const* tessellation,
		vector<vector<double> > const& tracers,double dm,
		Edge const& edge,int index,double dt,
		double time,ScalarInterpolation const* interp) const;

private:
  const Primitive outer_;
  RiemannSolver const& rs_;
  const vector<double> outer_tracer_;
};

#endif // FREEFLOW_HPP
