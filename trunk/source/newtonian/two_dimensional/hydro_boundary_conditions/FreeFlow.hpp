/*! \file FreeFlow.hpp
  \brief Free Flow Hydro Boundary Conditions
  \author Elad Steinberg
*/

#ifndef FREEFLOW_HPP
#define FREEFLOW_HPP 1

#include <string>
#include <iostream>
#include <fstream>
#include "../HydroBoundaryConditions.hpp"

//! \brief Free Flow Hydro Boundary Conditions
class FreeFlow: public HydroBoundaryConditions
{
public:
	/*!
	\brief Class constructor
	\param rs The Riemann solver
	\param mass_count Flag whether to keep track of the mass that was removed
	*/
  FreeFlow(RiemannSolver const& rs,bool mass_count=false);
	/*!
	\brief Class destructor
	*/
  ~FreeFlow(void);

  Conserved CalcFlux(Tessellation const* tessellation,
		     vector<Primitive> const& cells,Vector2D const& edge_velocity,
		     Edge const& edge,
		     SpatialReconstruction const* interp,double dt,
		     double time) const;

  Vector2D CalcEdgeVelocity(Tessellation const* tessellation,
			    vector<Vector2D> const& point_velocities,
			    Edge const& edge, double time) const;
  
  bool IsBoundary(Edge const& edge,Tessellation const* Data)const;
  
  bool IsGhostCell(int i,Tessellation const* Data) const;

  Primitive GetBoundaryPrimitive(Edge const& edge,
	  Tessellation const* Data,vector<Primitive> const& cells,
	  double time)const;

  vector<double> GetBoundaryTracers(Edge const& edge,
	  Tessellation const* Data,vector<vector<double> > const& tracers,double time)const;

  vector<double> CalcTracerFlux
  (Tessellation const* tessellation,vector<Primitive> const& cells,
   vector<vector<double> > const& tracers,double dm,
   Edge const& edge,int index,double dt,
   double time,
   SpatialReconstruction const* interp,Vector2D const& edge_velocity) const;

  /*!
  \brief Returns the ejected mass flux
  \return The mass flux
  */
  double GetMassFlux(void)const;
  /*!
  \brief Reads the ejected mass flux from a file
  \param loc The file's location
  */
  void ReadMassFlux(string loc)const;
  /*!
  \brief Writes the ejected mass flux to a file
  \param loc The file's location
  */
  void WriteMassFlux(string loc)const;
private:
  RiemannSolver const& rs_;
  bool mass_count_;
  mutable double mass_flux;
};

#endif // FREEFLOW_HPP
