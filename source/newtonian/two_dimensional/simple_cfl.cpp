#include "simple_cfl.hpp"
#include "hydrodynamics_2d.hpp"
#include <boost/foreach.hpp>
#ifdef RICH_MPI
#include <mpi.h>
#endif

SimpleCFL::SimpleCFL(const double cfl): cfl_(cfl) {}

namespace {
  class TimeStepCalculator: public LazyList<double>
  {
  public:

    TimeStepCalculator(const Tessellation& tess,
		       const vector<ComputationalCell>& cells,
		       const EquationOfState& eos,
		       const vector<Vector2D>& edge_velocities,
		TracerStickerNames const& tracerstickernames):
      tess_(tess), cells_(cells), 
      edge_velocities_(edge_velocities), eos_(eos), tracerstickernames_(tracerstickernames){}

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
	 cells_[i].tracers,
		tracerstickernames_.tracer_names);
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
	TracerStickerNames const& tracerstickernames_;
  };
}

double SimpleCFL::operator()(const Tessellation& tess,
			     const vector<ComputationalCell>& cells,
			     const EquationOfState& eos,
			     const vector<Vector2D>& point_velocities,
			     const double /*time*/,TracerStickerNames const& tracerstickernames) const
{
  double res =  cfl_*lazy_min(TimeStepCalculator(tess,cells,eos,point_velocities,tracerstickernames));
#ifdef RICH_MPI
  MPI_Allreduce(&res, &res, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
  return res;
}
