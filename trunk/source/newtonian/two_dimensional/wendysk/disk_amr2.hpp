/*! \file disk_amr2.hpp
  \brief AMR scheme for an accretion disk with spiral shocks
  \author Elad Steinberg
 */

#ifndef DISK_AMR_HPP
#define DISK_AMR_HPP 1

#include "../RemovalStrategy.hpp"
#include "../RefineStrategy.hpp"

//! \brief Coarsening scheme for an accretion disk with spiral shocks
class DiskRemove2: public RemovalStrategy
{
public:

  /*! \brief Class constructor
    \param inner_radius1 TBA
    \param inner_radius2 TBA
    \param outer_radius TBA
    \param total_specials Total number of non hydrodynamical cells
    \param s1 TBA
    \param s2 TBA
    \todo Add documentation
   */
  DiskRemove2(double inner_radius1,double inner_radius2,double outer_radius,
	  int total_specials,Vector2D s1,Vector2D s2);

  vector<int> CellsToRemove(Tessellation const* tess,
			    vector<Primitive> const& /*cells*/,
			    vector<vector<double> > const& /*tracers*/,
			    double /*time*/) const;

private:
  const double inner_radius1_;
  const double inner_radius2_;
  const double outer_radius_;
  const int total_specials_;
  const Vector2D s1_;
  const Vector2D s2_;
};

//! \brief Refinement scheme
class DiskRefine2: public RefineStrategy
{
public:

  /*! \brief Class constructor
    \param total_specials Total number of non hydrodynamical cells
    \param d_min Minimum density
    \param min_radius1 TBA
    \param min_radius2 TBA
    \param max_radius TBA
    \param max_volume Upper bound on cell volume
    \param s1 TBA
    \param s2 TBA
    \todo Add documentation
   */
  DiskRefine2(int total_specials,double d_min,double min_radius1,double min_radius2,
	  double max_radius,double max_volume,Vector2D s1,Vector2D s2);

  vector<int> CellsToRefine
  (Tessellation const* tess,
   vector<Primitive> const& cells,vector<vector<double> > const& /*tracers*/,
   double /*time*/,vector<Vector2D> & directions ,vector<int> const& Removed);

private:
  const int total_specials_;
  const double d_min_;
  const double min_radius1_;
  const double min_radius2_;
  const double max_radius_;
  const double max_volume_;
  const Vector2D s1_;
  const Vector2D s2_;
};

#endif // DISK_AMR_HPP
