#include "CourantFriedrichsLewy.hpp"
#include "../../misc/utils.hpp"
#include "../../misc/lazy_list.hpp"
#include <limits>
#ifdef RICH_MPI
#include <mpi.h>
#endif

CourantFriedrichsLewy::CourantFriedrichsLewy(double cfl, double SourceCFL, SourceTerm3D const& source, bool debug) :
	cfl_(cfl), sourcecfl_(SourceCFL), source_(source), debug_(debug), first_try_(true), dt_first_(-1)
{
	assert(cfl_ < 1 && "cfl number must be smaller than 1");
}

double CourantFriedrichsLewy::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
	const EquationOfState& eos, const vector<Vector3D>& face_velocities, const double /*time*/,
	TracerStickerNames const& tracerstickernames) const
{
	double res = 0.001*std::numeric_limits<double>::max();
	size_t N = tess.GetPointNo();
	size_t loc = 0;
	if (N > 0)
	{
		for (size_t i = 0; i < N; ++i)
		{
			double res_temp = 0;
			double c = 0;
#ifdef RICH_DEBUG
			try
			{
#endif
				c = eos.dp2c(cells[i].density, cells[i].pressure, cells[i].tracers,
					tracerstickernames.tracer_names);
#ifdef RICH_DEBUG
			}
			catch (UniversalError &eo)
			{
				eo.AddEntry("Error in CFL", 0);
				eo.AddEntry("Cell number", i);
				throw eo;
			}
#endif
			Vector3D const& v = cells[i].velocity;
			face_vec const& faces = tess.GetCellFaces(i);
			size_t Nloop = faces.size();
			for (size_t j = 0; j < Nloop; ++j)
			{
				Vector3D n = tess.Normal(faces[j]);
				n *= 1.0/fastabs(n);
				res_temp = fmax(res_temp, (c + std::abs(ScalarProd(n,v - face_velocities[faces[j]]))));
			}
			res_temp = tess.GetWidth(i) / res_temp;
			if (res_temp < res)
			{
				res = res_temp;
				loc = i;
			}
		}
	}
	res *= cfl_;
	res = 1.0 / std::max(source_.SuggestInverseTimeStep() / sourcecfl_, 1.0 / res);
	double old_res = res;
#ifdef RICH_MPI
	double new_res = 0;
	MPI_Allreduce(&res, &new_res, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	res = new_res;
#endif
	if (first_try_ && dt_first_ > 0)
	{
		res = dt_first_;
		first_try_ = false;
	}
	if (debug_)
	{
		if (1.000001*res > old_res)
		{
			Vector3D const& v = cells[loc].velocity;
			double c = eos.dp2c(cells[loc].density, cells[loc].pressure, cells[loc].tracers,
				tracerstickernames.tracer_names);
			std::cout << "Min dt, cell ID " << cells[loc].ID<<" width "<<tess.GetWidth(loc)<<" c "
				<<c<<" cell_loc "<<tess.GetMeshPoint(loc).x<<"," << tess.GetMeshPoint(loc).y << "," << tess.GetMeshPoint(loc).z
				<<" cell v "<<cells[loc].velocity.x<<"," << cells[loc].velocity.y << "," << cells[loc].velocity.z << std::endl;
			face_vec const& faces = tess.GetCellFaces(loc);
			size_t Nloop = faces.size();
			for (size_t j = 0; j < Nloop; ++j)
				std::cout << " face_vel " << fastabs(v - face_velocities[faces[j]])<<" ";
			std::cout<<std::endl;
		}
	}
	return res;
}

void CourantFriedrichsLewy::SetTimeStep(double dt)
{
	dt_first_ = dt;
	first_try_ = true;
}
