#ifndef DISK_AMR_HPP
#define DISK_AMR_HPP 1

#include "../RemovalStrategy.hpp"
#include "../RefineStrategy.hpp"

class DiskRemove: public RemovalStrategy
{
public:

  DiskRemove(double inner_radius,
	     double outer_radius,
	     int total_specials);

  vector<int> CellsToRemove(Tessellation const* tess,
			    vector<Primitive> const& /*cells*/,
			    vector<vector<double> > const& /*tracers*/,
			    double /*time*/) const;

private:
  const double inner_radius_;
  const double outer_radius_;
  const int total_specials_;
};

class DiskRefine: public RefineStrategy
{
public:

  DiskRefine(int total_specials,
	     double d_min,
	     double min_radius,double max_radius,
	     double max_volume);

  vector<int> CellsToRefine
  (Tessellation const* tess,
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
