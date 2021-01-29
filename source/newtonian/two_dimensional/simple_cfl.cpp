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
}

double SimpleCFL::operator()(const Tessellation& tess,
			     const vector<ComputationalCell>& cells,
			     const EquationOfState& eos,
			     const vector<Vector2D>& point_velocities,
			     const double /*time*/,TracerStickerNames const& tracerstickernames) const
{
  const vector<double> candidates =
    serial_generate<size_t, double>
    (create_range<size_t>(0, cells.size()),
     [&](size_t i)
     {return calc_local_time_step(tess,
				  cells,
				  point_velocities,
				  eos,
				  tracerstickernames,
				  i);});
  const double res =  cfl_*(*min_element(candidates.begin(),
					 candidates.end()));
#ifdef RICH_MPI
  double global_res=res;
  MPI_Allreduce(&res, &global_res, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  return global_res;
#else
  return res;
#endif
}
