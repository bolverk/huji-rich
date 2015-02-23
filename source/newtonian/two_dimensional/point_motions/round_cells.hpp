/*! \file round_cells.hpp
  \brief Adds small velocity correction to another point motion in order to make the cells rounder
  \author Elad Steinberg
 */

#ifndef ROUND_CELLS_HPP
#define ROUND_CELLS_HPP 1

#include "../point_motion.hpp"
#include "../hydrodynamics_2d.hpp"
#include "../OuterBoundary.hpp"
#include "../../../mpi/mpi_macro.hpp"

//! \brief Adds small velocity correction to another point motion in order to make the cells rounder
class RoundCells: public PointMotion
{
public:

  Vector2D CalcVelocity
  (int index, Tessellation const& tessellation,
   vector<Primitive> const& primitives,double /*time*/);

  /*! \brief Class constructor
    \param pm Base point motion
	\param hbc The hydro boundary conditions
    \param chi Chi (see equation 63 in arepo paper). Default value is 1 if coldflows is on.
    \param eta Eta (see equation 63 in arepo paper)
	\param coldflows The coldflows flag. Chooses if the fix velocity is set by the soundsspeed or by the time step
    \param innerNum The number of innerboundary cells
    \param outer pointer to outer boundary conditions, used to slow points so they don't overshoot the boundaries
  */
  RoundCells(PointMotion& pm,
	     HydroBoundaryConditions const& hbc,
	     double chi=0.15,
	     double eta=0.02,
	     bool coldflows=false,
	     int innerNum=0,
	     OuterBoundary const* outer=0);

  vector<Vector2D> calcAllVelocities
  (Tessellation const& tess,
   vector<Primitive> const& cells,
   double time,vector<CustomEvolution*> &cevolve,
   const vector<vector<double> >& tracers);

  /*!
  \brief Sets an external time step
  \param dt The time step
  */
  void SetExternalTimeStep(double dt);

private:
	void CorrectPointsOverShoot(vector<Vector2D> &v,double dt,
		Tessellation const& tess) const;

  PointMotion& pm_;
  HydroBoundaryConditions const& hbc_;
  const double chi_;
  const double eta_;
  const int inner_;
  OuterBoundary const* outer_;
  bool coldflows_;
  double lastdt_;
  bool evencall_;
  double external_dt_;

  RoundCells(const RoundCells& origin);
  RoundCells& operator=(const RoundCells& origin);
};

#endif // ROUND_CELLS_HPP
