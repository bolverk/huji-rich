#include "simple_cfl.hpp"
#include "hydrodynamics_2d.hpp"

SimpleCFL::SimpleCFL(const double cfl): cfl_(cfl) {}

namespace {
  class TimeStepCalculator: public LazyList<double>
  {
  public:

    TimeStepCalculator(const Tessellation& tess,
		       const vector<ComputationalCell>& cells,
		       const EquationOfState& eos,
		       const vector<Vector2D>& point_velocities):
      tess_(tess), cells_(cells), 
      point_velocities_(point_velocities), eos_(eos) {}

    size_t size(void) const 
    {
      return cells_.size();
    }

    double operator[](size_t i) const
    {
      return tess_.GetWidth(static_cast<int>(i))/
	(eos_.dp2c(cells_[i].density, cells_[i].pressure, cells_[i].tracers)+
	 abs(cells_[i].velocity)+
	 abs(point_velocities_[i]));
    }

  private:
    const Tessellation& tess_;
    const vector<ComputationalCell>& cells_;
    const vector<Vector2D>& point_velocities_;
    const EquationOfState& eos_;
  };
}

double SimpleCFL::operator()(const Tessellation& tess,
			     const vector<ComputationalCell>& cells,
			     const EquationOfState& eos,
			     const vector<Vector2D>& point_velocities,
			     const double /*time*/) const
{
  return cfl_*lazy_min(TimeStepCalculator(tess,cells,eos,point_velocities));
}
