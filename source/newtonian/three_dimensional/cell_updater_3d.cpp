#include "cell_updater_3d.hpp"
#include <iostream>

CellUpdater3D::~CellUpdater3D(void) {}

namespace
{
	class SolveVelocity
	{

	public:
		double a0_, a1_, a2_, a3_, a4_;

		SolveVelocity(double a0, double a1, double a2, double a3, double a4) :a0_(a0), a1_(a1), a2_(a2), a3_(a3), a4_(a4) {}

		double operator()(const double v) const
		{
			return v * v*(a4_*v*v + a3_ * v + a2_) + v * a1_ + a0_;
			//return a4_ * v*v*v*v + a3_ * v*v*v + a2_ * v*v + a1_ * v + a0_;
		}

		double Deriv(const double v) const
		{
			return v * v*(4 * a4_*v + a3_ * 3) + 2 * a2_*v + a1_;
			//return 4 * a4_*v*v*v + 3 * a3_*v*v + 2 * a2_ * v + a1_;
		}
	};

	double DoNewtonRapshon(SolveVelocity const& solve, double val)
	{
		size_t counter = 1;
		double f0 = solve(val);
		double new_val = val - f0 / solve.Deriv(val);
		while (std::abs(new_val - val) > 1e-12 && (std::abs(f0) > (1e-12*solve.a0_)))
		{
			++counter;
			val = new_val;
			f0 = solve(val);
			new_val = std::min(1.0, val - f0 / solve.Deriv(val));
			if (counter > 99)
			{
				throw UniversalError("Bad convergence in simple cell updater, too mant iterations in finding velocity");
			}
		}
		return new_val;
	}
}

double GetVelocity(Conserved3D const& cell, double G)
{
	double M = std::sqrt(ScalarProd(cell.momentum, cell.momentum));
	// Add rest mass energy
	double E = cell.energy + cell.mass;
	SolveVelocity tosolve(M*M, -2 * G*M*E, G*G*E*E + 2 * (G - 1)*M*M - (G - 1)*(G - 1)*cell.mass*cell.mass, -2 * G*(G - 1)*M*E, (G - 1)*(G - 1)*(cell.mass*cell.mass + M * M));

	double vmin = (1e6*M < cell.mass) ? 0 : (G*E - std::sqrt((G*E)*(G*E) - 4 * (G - 1)*M*M)) / (2 * M*(G - 1));
	vmin = std::min(vmin, 0.999);
	if ((G*E)*(G*E) - 4 * (G - 1)*M*M < 0)
		vmin = 0;
	double vmax = std::min(1.0, M / E + 1e-6);
	double res = 0;
	try
	{
		res = DoNewtonRapshon(tosolve, 0.5*(vmin + vmax));
	}
	catch (UniversalError &eo)
	{
		eo.AddEntry("Mass", cell.mass);
		eo.AddEntry("Mx", cell.momentum.x);
		eo.AddEntry("My", cell.momentum.y);
		eo.AddEntry("Mz", cell.momentum.z);
		eo.AddEntry("Energy", cell.energy);
		eo.AddEntry("Enthalpy", cell.internal_energy);
		DisplayError(eo);
		throw eo;
	}
	return res;
}