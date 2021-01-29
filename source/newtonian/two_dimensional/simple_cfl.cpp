#include "simple_cfl.hpp"
#include "hydrodynamics_2d.hpp"
#include <boost/foreach.hpp>
#include "../../misc/serial_generate.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif

SimpleCFL::SimpleCFL(const double cfl): cfl_(cfl) {}

namespace {

  double calc_local_time_step
  (const Tessellation& tess,
   const vector<ComputationalCell>&  cells,
   const vector<Vector2D>& edge_velocities,
   const EquationOfState& eos,
   const TracerStickerNames& tsn,
   size_t i)
  {
    const double radius = tess.GetWidth(static_cast<int>(i));
    const ComputationalCell& cell = cells.at(i);
    const double c = eos.dp2c(cell.density,
			      cell.pressure,
			      cell.tracers,
			      tsn.tracer_names);
    const Vector2D v = cell.velocity;
    const vector<double> candidates =
      serial_generate<int, double>
      (tess.GetCellEdges(static_cast<int>(i)),
       [&](int j)
       {
	 const Vector2D ve = edge_velocities.at(static_cast<size_t>(j));
	 const Edge& edge = tess.GetEdge(j);
	 const Vector2D n =
	   normalize(tess.GetMeshPoint(edge.neighbors.second)-
		     tess.GetMeshPoint(edge.neighbors.first));
	 return radius/(c+abs(ScalarProd(n,v-ve)));
       });
    return *min_element(candidates.begin(),
			candidates.end());
  }

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
      /*
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
	Edge const& edge = tess_.GetEdge(index);
	const Vector2D n =
	  normalize(tess_.GetMeshPoint(edge.neighbors.second) - tess_.GetMeshPoint(edge.neighbors.first));
	res = fmax
	  (res,
	   (c+std::abs(ScalarProd(n,(v-ve))))/radius);
      }
      return 1.0/res;
      */
      return calc_local_time_step(tess_,
				  cells_,
				  edge_velocities_,
				  eos_,
				  tracerstickernames_,
				  i);
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
  double global_res=res;
  MPI_Allreduce(&res, &global_res, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  return global_res;
#else
  return res;
#endif
}
