/*! \file disk_amr.hpp
  \brief AMR scheme for accretion disk with spiral shocks
  \author Almog Yalinewich
 */

#ifndef DISK_AMR_HPP
#define DISK_AMR_HPP 1

#include "../RemovalStrategy.hpp"
#include "../RefineStrategy.hpp"
#include <cmath>

//! \brief Coarsening scheme
class DiskRemove: public RemovalStrategy
{
public:

  /*! \brief 
    \param inner_radius Inner bound of computational domain
    \param outer_radius Outer bound of computational domain
    \param total_specials Total number of non hydrodynamical cells
    \param MinVol Lower bound on cell volume
   */
  DiskRemove(double inner_radius,
	     double outer_radius,
	     int total_specials,double MinVol);

  vector<int> CellsToRemove(Tessellation const& tess,
			    vector<Primitive> const& /*cells*/,
			    vector<vector<double> > const& /*tracers*/,
			    double /*time*/) const;

private:
  const double inner_radius_;
  const double outer_radius_;
  const int total_specials_;
  const double MinVol_;
};

//! \brief Coarsening scheme for accretion disk with spiral shocks
class DiskRefine: public RefineStrategy
{
public:

  /*! \brief Class constructor
    \param total_specials Total number of non hydrodynamical cells
    \param d_min Minimum density
    \param min_radius Inner boundary of the domain
    \param max_radius Outer boundary of the domain
    \param max_volume Upper bound on cell volume
   */
  DiskRefine(int total_specials,
	     double d_min,
	     double min_radius,double max_radius,
	     double max_volume);

  vector<int> CellsToRefine
  (Tessellation const& tess,
   vector<Primitive> const& cells,vector<vector<double> > const& /*tracers*/,
   double /*time*/,vector<Vector2D> & directions ,vector<int> const& Removed);

private:
  const int total_specials_;
  const double d_min_;
  const double min_radius_;
  const double max_radius_;
  const double max_volume_;
};

#endif // DISK_AMR_HPP
