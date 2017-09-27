#include "CourantFriedrichsLewy.hpp"
#include "../../misc/utils.hpp"
#include "../../misc/lazy_list.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif

CourantFriedrichsLewy::CourantFriedrichsLewy(double cfl, SourceTerm3D const& source):
cfl_(cfl),source_(source),first_try_(true),dt_first_(-1)
{
	assert(cfl_<1 && "cfl number must be smaller than 1");
}

namespace 
{
	class TimeStepBoundCalculator: public LazyList<double>
	{
	public:

		TimeStepBoundCalculator(const Tessellation3D& tess,const vector<ComputationalCell3D>& cells,
			const EquationOfState& eos, const vector<Vector3D>& face_velocities,TracerStickerNames const&
			tracerstickernames):
		tess_(tess), cells_(cells), eos_(eos),face_velocities_(face_velocities),tracerstickernames_(tracerstickernames) {}

	  size_t size(void) const
		{
			return tess_.GetPointNo();
		}

		double operator[](size_t i) const
		{
			double res = 0;
			const double c = eos_.dp2c(cells_[i].density,cells_[i].pressure,cells_[i].tracers,
				tracerstickernames_.tracer_names);
			const Vector3D v = cells_.at(i).velocity;
			vector<size_t> const& faces = tess_.GetCellFaces(i);
			size_t Nloop = faces.size();
			for (size_t j = 0; j < Nloop; ++j)
				res = fmax(res,(c + abs(v - face_velocities_[faces[j]])));
			return tess_.GetWidth(i) / res;
		}

	private:
		const Tessellation3D& tess_;
		const vector<ComputationalCell3D>& cells_;
		const EquationOfState& eos_;
		const vector<Vector3D>& face_velocities_;
		TracerStickerNames const& tracerstickernames_;
	};
}

double CourantFriedrichsLewy::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells, 
	const EquationOfState& eos,const vector<Vector3D>& face_velocities, const double /*time*/,
	TracerStickerNames const& tracerstickernames) const
{
	double res = cfl_*lazy_min(TimeStepBoundCalculator(tess,cells,eos,face_velocities,tracerstickernames));
	res = 1.0 / std::max(source_.SuggestInverseTimeStep()/cfl_, 1.0 / res);
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
	return res;
}

void CourantFriedrichsLewy::SetTimeStep(double dt)
{
	dt_first_ = dt;
}