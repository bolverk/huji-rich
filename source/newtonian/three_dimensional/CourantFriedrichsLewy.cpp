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

		TimeStepBoundCalculator(const Tessellation3D& tess,
			const vector<ComputationalCell>& cells,
			const EquationOfState& eos):
		tess_(tess), cells_(cells), eos_(eos) {}

	  size_t size(void) const
		{
			return cells_.size();
		}

		double operator[](size_t i) const
		{
			return tess_.GetWidth(i)/
				(abs(cells_[i].velocity)+
				eos_.dp2c(cells_[i].density,cells_[i].pressure));
		}

	private:
		const Tessellation3D& tess_;
		const vector<ComputationalCell>& cells_;
		const EquationOfState& eos_;
	};
}

double CourantFriedrichsLewy::operator()
	(const Tessellation3D& tess,
	const vector<ComputationalCell>& cells,
	const EquationOfState& eos) const
{
	return cfl_*lazy_min(TimeStepBoundCalculator(tess,cells,eos));
}
