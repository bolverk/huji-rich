#ifndef ROUND_CELLS_HPP
#define ROUND_CELLS_HPP 1

#include "../point_motion.hpp"
#include "../hydrodynamics_2d.hpp"
#include "../OuterBoundary.hpp"

//! \brief Adds small velocity correction to another point motion in order to make the cells rounder
class RoundCells: public PointMotion
{
public:
	
  Vector2D CalcVelocity
  (int index, Tessellation const* tessellation,
   vector<Primitive> const& primitives,double /*time*/);
	
  /*! \brief Class constructor
    \param pm Base point motino
    \param chi Chi (see equation 63 in arepo paper)
    \param eta Eta (see equation 63 in arepo paper)
    \param innerNum The number of innerboundary cells
    \param outer pointer to outer boundary conditions, used to slow points so they don't overshoot the boundaries
  */
  RoundCells(PointMotion& pm, double chi=1.0,
	     double eta=0.1,int innerNum=0,OuterBoundary const* outer=0);

  /*! \brief Calculates the velocities of all cells
    \param tess Tessellation
    \param cells Hydrodynamic cells
    \param time Time
   */
  vector<Vector2D> calcAllVelocities
  (Tessellation const* tess,
   vector<Primitive> const& cells,
   double time);

  /*! \brief Sets the option for cold flows dependence
    \param hbc The hydro boundary conditions
  */
  void SetColdFlows(HydroBoundaryConditions *hbc);

  /*!
  \brief Sets an external time step
  \param dt The time step
  */
  void SetExternalTimeStep(double dt);

private:
	void CorrectPointsOverShoot(vector<Vector2D> &v,double dt,
		Tessellation const* tess) const;

  PointMotion& pm_;
  const double chi_;
  const double eta_;
  const int inner_;
  OuterBoundary const* outer_;
  bool coldflows_;
  HydroBoundaryConditions *hbc_;
  double lastdt_;
  bool evencall_;
  double external_dt_;

  RoundCells(const RoundCells& origin);
  RoundCells& operator=(const RoundCells& origin);
};

#endif // ROUND_CELLS_HPP
