#include "CourantFriedrichsLewy.hpp"
#include "../../misc/utils.hpp"
#include "../../misc/lazy_list.hpp"

CourantFriedrichsLewy::CourantFriedrichsLewy(double cfl):
cfl_(cfl)
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
	return cfl_*lazy_min(TimeStepBoundCalculator(tess,cells,eos,face_velocities,tracerstickernames));
}
