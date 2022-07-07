#include "CourantFriedrichsLewy.hpp"
#include "../../misc/utils.hpp"
#include "../../misc/lazy_list.hpp"
#include <limits>
#ifdef RICH_MPI
#include <mpi.h>
#endif

CourantFriedrichsLewy::CourantFriedrichsLewy(double cfl, double SourceCFL, SourceTerm3D const& source, 
	std::vector<std::string> no_calc, bool debug) :
	cfl_(cfl), sourcecfl_(SourceCFL), source_(source), no_calc_(no_calc), debug_(debug), first_try_(true), dt_first_(-1), last_time_(-10000)
{
	assert(cfl_ < 1 && "cfl number must be smaller than 1");
}

double CourantFriedrichsLewy::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
	const EquationOfState& eos, const vector<Vector3D>& face_velocities, const double time) const
{
	double res = 0.001 * std::numeric_limits<double>::max();
	size_t const N = tess.GetPointNo();
	size_t loc = 0;
	size_t const N_no_calc = no_calc_.size();
	std::vector<size_t> no_calc_indeces(N_no_calc);
	for(size_t i = 0; i <N_no_calc; ++i)
		no_calc_indeces[i] = binary_index_find(ComputationalCell3D::stickerNames, no_calc_[i]);
	if (N > 0)
	{
		for (size_t i = 0; i < N; ++i)
		{
			bool cont = false;
			for(size_t j = 0; j < N_no_calc; ++j)
				if(cells[i].stickers[no_calc_indeces[j]])
					cont = true;
			if(cont)
				continue;
			double res_temp = 0;
			double c = 0;
#ifdef RICH_DEBUG
			try
			{
#endif
				c = eos.dp2c(cells[i].density, cells[i].pressure, cells[i].tracers,
					     ComputationalCell3D::tracerNames);
#ifdef RICH_DEBUG
			}
			catch (UniversalError& eo)
			{
				eo.addEntry("Error in CFL", 0);
				eo.addEntry("Cell number", i);
				throw eo;
			}
#endif
			Vector3D const& v = cells[i].velocity;
			face_vec const& faces = tess.GetCellFaces(i);
			size_t const Nloop = faces.size();
			double max_face_area = 0;
			for (size_t j = 0; j < Nloop; ++j)
			{
				Vector3D n = tess.Normal(faces[j]);
				n *= 1.0 / fastabs(n);
				res_temp = fmax(res_temp, (c + std::abs(ScalarProd(n, v - face_velocities[faces[j]]))));
				max_face_area = std::max(max_face_area, tess.GetArea(faces[j]));
			}
			double cell_effective_radius = std::min(tess.GetWidth(i), tess.GetVolume(i) / 
				max_face_area);
			res_temp = cell_effective_radius / res_temp;
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
	if (first_try_ && dt_first_ > 0 || (last_time_ == time && dt_first_ > 0))
	{
		res = std::min(res, dt_first_);
		first_try_ = false;
		if (close2zero(last_time_ - time))
			dt_first_ = -1;
	}
	if (debug_)
	{
		if (1.0000001 * res > old_res)
		{
			Vector3D const& v = cells[loc].velocity;
			double c = eos.dp2c(cells[loc].density, cells[loc].pressure, cells[loc].tracers, ComputationalCell3D::tracerNames);
			std::cout << "Min dt="<<res<<", cell ID " << cells[loc].ID << " width " << tess.GetWidth(loc) << " c "
				<< c << " cell v " << cells[loc].velocity.x << "," << cells[loc].velocity.y << "," << cells[loc].velocity.z << std::endl;
			face_vec const& faces = tess.GetCellFaces(loc);
			size_t Nloop = faces.size();
			for (size_t j = 0; j < Nloop; ++j)
			{
				std::cout << " face_vel " << fastabs(v - face_velocities[faces[j]]) << " ";
				Vector3D p1 = tess.GetMeshPoint(tess.GetFaceNeighbors(faces[j]).first);
				Vector3D p2 = tess.GetMeshPoint(tess.GetFaceNeighbors(faces[j]).second);
				std::cout<<"p1="<<p1.x<<","<<p1.y<<","<<p1.z<<" p2="<<p2.x<<","<<p2.y<<","<<p2.z<<std::endl;
				p1 = cells[tess.GetFaceNeighbors(faces[j]).first].velocity;
				p2 = cells[tess.GetFaceNeighbors(faces[j]).second].velocity;
				std::cout<<"v1="<<p1.x<<","<<p1.y<<","<<p1.z<<" v2="<<p2.x<<","<<p2.y<<","<<p2.z<<
				" A "<<tess.GetArea(faces[j])<<std::endl;
			}
		}
	}
	last_time_ = time;
	return res;
}

void CourantFriedrichsLewy::SetTimeStep(double dt)
{
	dt_first_ = dt;
	first_try_ = true;
}
