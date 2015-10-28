#include "simple_cfl.hpp"
#include "hydrodynamics_2d.hpp"
#include <boost/foreach.hpp>
#ifdef RICH_MPI
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#endif

SimpleCFL::SimpleCFL(const double cfl): cfl_(cfl) {}

namespace {
  class TimeStepCalculator: public LazyList<double>
  {
  public:

    TimeStepCalculator(const Tessellation& tess,
		       const vector<ComputationalCell>& cells,
		       const EquationOfState& eos,
		       const vector<Vector2D>& edge_velocities):
      tess_(tess), cells_(cells), 
      edge_velocities_(edge_velocities), eos_(eos) {}

    size_t size(void) const 
    {
		return static_cast<size_t>(tess_.GetPointNo());
    }

    double operator[](size_t i) const
    {
      double res = 0;
      const double radius = tess_.GetWidth(static_cast<int>(i));
      const double c = eos_.dp2c
	(cells_[i].density,
	 cells_[i].pressure,
	 cells_[i].tracers);
      const Vector2D v = cells_.at(i).velocity;
      BOOST_FOREACH
	(int index,
	 tess_.GetCellEdges(static_cast<int>(i))){
	const Vector2D ve = edge_velocities_.at(static_cast<size_t>(index));
	res = fmax
	  (res,
	   (c+abs(v-ve))/radius);
      }
      return 1.0/res;
    }

  private:
    const Tessellation& tess_;
    const vector<ComputationalCell>& cells_;
    const vector<Vector2D>& edge_velocities_;
    const EquationOfState& eos_;
  };
}

double SimpleCFL::operator()(const Tessellation& tess,
			     const vector<ComputationalCell>& cells,
			     const EquationOfState& eos,
			     const vector<Vector2D>& point_velocities,
			     const double /*time*/) const
{
  double res =  cfl_*lazy_min(TimeStepCalculator(tess,cells,eos,point_velocities));
#ifdef RICH_MPI
  const boost::mpi::communicator world;
  res = boost::mpi::all_reduce(world, res, boost::mpi::minimum<double>());
#endif
  return res;
}
