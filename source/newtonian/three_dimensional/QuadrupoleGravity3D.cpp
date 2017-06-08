#include "QuadrupoleGravity3D.hpp"
#include <iostream>
#include "../../misc/simple_io.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif

namespace
{
	double LinearInterp(const vector<double> &x, const vector<double> &y, double xi, size_t &index)
	{
		vector<double>::const_iterator it = std::upper_bound(x.begin(), x.end(), xi);
		assert(it != x.end() && it != x.begin() &&
			"X out of range in Linear Interp");
		index = std::max(static_cast<size_t>(it - x.begin()),static_cast<size_t>(1));
		if (*it == xi)
			return y[static_cast<std::size_t>(it - x.begin())];

		return y[static_cast<std::size_t>(it - x.begin())] + (xi - *it)*
			(y[static_cast<std::size_t>(it - 1 - x.begin())] -
				y[static_cast<std::size_t>(it - x.begin())]) / (*(it - 1) - *it);
	}

	double LinearInterp2(const vector<double> &x, const vector<double> &y, double xi, size_t index)
	{
		double x0 = x[index];
		return y[index] + (xi - x0)*	(y[index - 1] - y[index]) / (x[index - 1] - x0);
	}

	Vector3D CalcQuadGrad(vector<double> const& r_list, vector<double> const& Q20In, vector<double> const& Q21InR,
		vector<double> const& Q22InR, vector<double> const& Q20Out, vector<double> const& Q21OutR, vector<double> const& Q22OutR,
		vector<double> const& Q21InI, vector<double> const& Q22InI, vector<double> const& Q21OutI, vector<double> const& Q22OutI,
		Vector3D point)
	{
		Vector3D acc;
		size_t index = 0;
		double r = abs(point);
		double Phi = LinearInterp(r_list, Q20In, r, index);
		double small = 1e-5;
		double dx = small*(r_list[index] - r_list.at(index-1));
		point.x += 0.5*dx;
		r = abs(point);
		Phi = -LinearInterp(r_list, Q20In, r, index)*(2 * point.z*point.z - point.x*point.x - point.y*point.y) / (r*r*r*r*r);
		point.x -= dx;
		r = abs(point);
		Phi += LinearInterp2(r_list, Q20In, r, index)*(2 * point.z*point.z - point.x *point.x - point.y * point.y) / (r*r*r*r*r);
		point.x += dx;
		r = abs(point);
		Phi -= LinearInterp2(r_list, Q20Out, r, index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y);
		point.x -= dx;
		r = abs(point);
		Phi += LinearInterp2(r_list, Q20Out, r, index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y);
		point.x += dx;
		r = abs(point);
		Phi += 2 * LinearInterp2(r_list, Q21InR, r, index)*point.x*point.z / (r*r*r*r*r);
		point.x -= dx;
		r = abs(point);
		Phi -= 2 * LinearInterp2(r_list, Q21InR, r, index)*point.x*point.z / (r*r*r*r*r);
		point.x += dx;
		r = abs(point);
		Phi += 2 * LinearInterp2(r_list, Q21InI, r, index)*point.y*point.z / (r*r*r*r*r);
		point.x -= dx;
		r = abs(point);
		Phi -= 2 * LinearInterp2(r_list, Q21InI, r, index)*point.y*point.z / (r*r*r*r*r);
		point.x += dx;
		r = abs(point);
		Phi += 2 * LinearInterp2(r_list, Q21OutR, r, index)*point.x*point.z;
		point.x -= dx;
		r = abs(point);
		Phi -= 2 * LinearInterp2(r_list, Q21OutR, r, index)*point.x*point.z;
		point.x += dx;
		r = abs(point);
		Phi += 2 * LinearInterp2(r_list, Q21OutI, r, index)*point.y*point.z;
		point.x -= dx;
		r = abs(point);
		Phi -= 2 * LinearInterp2(r_list, Q21OutI, r, index)*point.y*point.z;
		point.x += dx;
		r = abs(point);
		Phi -= 2 * LinearInterp2(r_list, Q22InR, r, index)*(point.x*point.x - point.y*point.y) / (r*r*r*r*r);
		point.x -= dx;
		r = abs(point);
		Phi += 2 * LinearInterp2(r_list, Q22InR, r, index)*(point.x*point.x - point.y*point.y) / (r*r*r*r*r);
		point.x += dx;
		r = abs(point);
		Phi -= 4 * LinearInterp2(r_list, Q22InI, r, index)*point.y*point.x / (r*r*r*r*r);
		point.x -= dx;
		r = abs(point);
		Phi += 4 * LinearInterp2(r_list, Q22InI, r, index)*point.y*point.x / (r*r*r*r*r);
		point.x += dx;
		r = abs(point);
		Phi -= 2 * LinearInterp2(r_list, Q22OutR, r, index)*(point.x*point.x - point.y*point.y);
		point.x -= dx;
		r = abs(point);
		Phi += 2 * LinearInterp2(r_list, Q22OutR, r, index)*(point.x*point.x - point.y*point.y);
		point.x += dx;
		r = abs(point);
		Phi -= 4 * LinearInterp2(r_list, Q22OutI, r, index)*point.y*point.x;
		point.x -= dx;
		r = abs(point);
		Phi += 4 * LinearInterp2(r_list, Q22OutI, r, index)*point.y*point.x;

		acc.x = Phi / dx;
		point.x += 0.5*dx;

		point.y += 0.5*dx;
		r = abs(point);
		Phi = -LinearInterp2(r_list, Q20In, r, index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y) / (r*r*r*r*r);
		Phi -= LinearInterp2(r_list, Q20Out, r, index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y);
		Phi += 2 * LinearInterp2(r_list, Q21InR, r, index)*point.x*point.z / (r*r*r*r*r);
		Phi += 2 * LinearInterp2(r_list, Q21InI, r, index)*point.y*point.z / (r*r*r*r*r);
		Phi += 2 * LinearInterp2(r_list, Q21OutR, r, index)*point.x*point.z;
		Phi += 2 * LinearInterp2(r_list, Q21OutI, r, index)*point.y*point.z;
		Phi -= 2 * LinearInterp2(r_list, Q22InR, r, index)*(point.x*point.x - point.y*point.y) / (r*r*r*r*r);
		Phi -= 4 * LinearInterp2(r_list, Q22InI, r, index)*point.y*point.x / (r*r*r*r*r);
		Phi -= 2 * LinearInterp2(r_list, Q22OutR, r, index)*(point.x*point.x - point.y*point.y);
		Phi -= 4 * LinearInterp2(r_list, Q22OutI, r, index)*point.y*point.x;

		point.y -= dx;
		r = abs(point);

		Phi += LinearInterp2(r_list, Q20In, r, index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y) / (r*r*r*r*r);
		Phi += LinearInterp2(r_list, Q20Out, r, index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y);
		Phi -= 2 * LinearInterp2(r_list, Q21InR, r, index)*point.x*point.z / (r*r*r*r*r);
		Phi -= 2 * LinearInterp2(r_list, Q21InI, r, index)*point.y*point.z / (r*r*r*r*r);
		Phi -= 2 * LinearInterp2(r_list, Q21OutR, r, index)*point.x*point.z;
		Phi -= 2 * LinearInterp2(r_list, Q21OutI, r, index)*point.y*point.z;
		Phi += 2 * LinearInterp2(r_list, Q22InR, r, index)*(point.x*point.x - point.y*point.y) / (r*r*r*r*r);
		Phi += 4 * LinearInterp2(r_list, Q22InI, r, index)*point.y*point.x / (r*r*r*r*r);
		Phi += 2 * LinearInterp2(r_list, Q22OutR, r, index)*(point.x*point.x - point.y*point.y);
		Phi += 4 * LinearInterp2(r_list, Q22OutI, r, index)*point.y*point.x;
		acc.y = Phi / dx;
		point.y += 0.5*dx;

		point.z += 0.5*dx;
		r = abs(point);
		Phi = -LinearInterp2(r_list, Q20In, r, index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y) / (r*r*r*r*r);
		Phi -= LinearInterp2(r_list, Q20Out, r, index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y);
		Phi += 2 * LinearInterp2(r_list, Q21InR, r, index)*point.x*point.z / (r*r*r*r*r);
		Phi += 2 * LinearInterp2(r_list, Q21InI, r, index)*point.y*point.z / (r*r*r*r*r);
		Phi += 2 * LinearInterp2(r_list, Q21OutR, r, index)*point.x*point.z;
		Phi += 2 * LinearInterp2(r_list, Q21OutI, r, index)*point.y*point.z;
		Phi -= 2 * LinearInterp2(r_list, Q22InR, r, index)*(point.x*point.x - point.y*point.y) / (r*r*r*r*r);
		Phi -= 4 * LinearInterp2(r_list, Q22InI, r, index)*point.y*point.x / (r*r*r*r*r);
		Phi -= 2 * LinearInterp2(r_list, Q22OutR, r, index)*(point.x*point.x - point.y*point.y);
		Phi -= 4 * LinearInterp2(r_list, Q22OutI, r, index)*point.y*point.x;

		point.z -= dx;
		r = abs(point);

		Phi += LinearInterp2(r_list, Q20In, r, index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y) / (r*r*r*r*r);
		Phi += LinearInterp2(r_list, Q20Out, r, index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y);
		Phi -= 2 * LinearInterp2(r_list, Q21InR, r, index)*point.x*point.z / (r*r*r*r*r);
		Phi -= 2 * LinearInterp2(r_list, Q21InI, r, index)*point.y*point.z / (r*r*r*r*r);
		Phi -= 2 * LinearInterp2(r_list, Q21OutR, r, index)*point.x*point.z;
		Phi -= 2 * LinearInterp2(r_list, Q21OutI, r, index)*point.y*point.z;
		Phi += 2 * LinearInterp2(r_list, Q22InR, r, index)*(point.x*point.x - point.y*point.y) / (r*r*r*r*r);
		Phi += 4 * LinearInterp2(r_list, Q22InI, r, index)*point.y*point.x / (r*r*r*r*r);
		Phi += 2 * LinearInterp2(r_list, Q22OutR, r, index)*(point.x*point.x - point.y*point.y);
		Phi += 4 * LinearInterp2(r_list, Q22OutI, r, index)*point.y*point.x;
		acc.z = Phi / dx;

		return acc;
	}
}

QuadrupoleGravity3D::QuadrupoleGravity3D(size_t res, double smoothlength,bool output) :edges_(vector<double>(res+1,-1)),
smoothlength_(smoothlength),output_(output),potential(vector<double>()){}

void QuadrupoleGravity3D::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
	const vector<Conserved3D>& /*fluxes*/, const double /*time*/, TracerStickerNames const& /*tracerstickernames*/,
	vector<Vector3D> &acc) const
{
	size_t N = tess.GetPointNo();
	Vector3D CM;
	double M = 0;
	for (size_t i = 0; i < N; ++i)
	{
		double m = tess.GetVolume(i)*cells[i].density;
		M += m;
		CM += tess.GetCellCM(i)*m;
	}
	double Mtot = M;
	Vector3D CMtot = CM;
#ifdef RICH_MPI
	MPI_Allreduce(&M, &Mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&CM.x, &CMtot.x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&CM.y, &CMtot.y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&CM.z, &CMtot.z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	CMtot *= (1.0 / Mtot);

	// GetCM, assume center for now
	size_t Nmid = edges_.size()-1;
	vector<double> mr2(Nmid+1, 0), Mrtot(Nmid + 1, 0), dr_all(Nmid-1, 0);
	if (edges_[0] < 0)
	{
		edges_[0] = 0;
		edges_.back() = abs(tess.GetBoxCoordinates().first - tess.GetBoxCoordinates().second);
		for (size_t i = 1; i < Nmid; ++i)
			edges_[i] = static_cast<double>(i)*edges_.back() / static_cast<double>(Nmid);
	}
	for (size_t j = 0; j < 3; ++j)
	{
		Mrtot.assign(Nmid + 1, 0);
		mr2.assign(Nmid + 1, 0);
		for (size_t i = 0; i < N; ++i)
		{
			double r_temp = abs(tess.GetMeshPoint(i) - CM);
			size_t index = static_cast<size_t>(std::lower_bound(edges_.begin(), edges_.end(), r_temp) - edges_.begin());
			mr2[index] += tess.GetVolume(i)*cells[i].density/(r_temp*r_temp+smoothlength_*smoothlength_);
		}
		vector<double> temp(Nmid + 1, 0);
#ifdef RICH_MPI
		MPI_Allreduce(&mr2[0], &temp[0], static_cast<int>(Nmid+1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		mr2 = temp;
#endif
		for (size_t i = 1; i < Nmid + 1; ++i)
		{
			Mrtot[i] = Mrtot[i - 1] + mr2[i];
			if (Mrtot[i] < (1 + 1e-7)*Mrtot[i - 1])
				Mrtot[i] = (1 + 1e-7)*Mrtot[i - 1];
		}
		double dr = Mrtot[Nmid]/static_cast<double>(Nmid);
		for (size_t i = 0; i < (Nmid-1); ++i)
			dr_all[i] = static_cast<double>(i + 1)*dr;
		temp = edges_;
		for (size_t i = 1; i < Nmid; ++i)
			edges_[i] = LinearInterpolation(Mrtot, temp, dr_all[i-1]);
	}

	size_t resolution = Nmid;
	vector<double> m_radius(resolution+1,0), QIn20(resolution+1,0), QIn21_real(resolution+1,0),
		QIn21_I(resolution+1,0), QOut20(resolution+1,0), QOut21_real(resolution+1,0), QOut21_I(resolution+1,0),
		QOut22_I(resolution+1,0), QOut22_real(resolution+1,0), QIn22_I(resolution+1,0), QIn22_real(resolution+1,0);
	double A21 = 3.0 / 2.0, A22 = 0.25*A21;
	// Find mass and Q moments in radius
	Vector3D diff;

	for (size_t i = 0; i < N; ++i)
	{
		diff = tess.GetMeshPoint(i);
		diff -= CMtot;
		double m = tess.GetVolume(i)*cells[i].density;
		double r = abs(diff);
		size_t whole = static_cast<size_t>(std::lower_bound(edges_.begin(), edges_.end(), r) - edges_.begin());
		//double fraction = 1 - (edges_[whole]-r) / (edges_[whole]-edges_[whole-1]);
		double fraction = 1;
		double r_inv = 1.0 / (r*r*r*r*r+smoothlength_*smoothlength_*smoothlength_*smoothlength_*smoothlength_);
		m_radius[whole] += m*fraction;
		QIn20[whole] += 0.25*m*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*fraction;
		QIn21_I[whole] -= m*A21*diff.y*diff.z*fraction;
		QIn21_real[whole] -= m*A21*diff.z*diff.x*fraction;
		QIn22_real[whole] += m*A22*(diff.x*diff.x - diff.y*diff.y)*fraction;
		QIn22_I[whole] += 2 * m*A22*(diff.x*diff.y)*fraction;
		whole--;
		QOut20[whole] += m*0.25*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*r_inv*fraction;
		QOut21_I[whole] -= m*A21*diff.y*diff.z*r_inv*fraction;
		QOut21_real[whole] -= m*A21*diff.z*diff.x*r_inv*fraction;
		QOut22_real[whole] += m*A22*(diff.x*diff.x - diff.y*diff.y)*r_inv*fraction;
		QOut22_I[whole] += 2 * m*A22*(diff.x*diff.y)*r_inv*fraction;
		/*m_radius[whole - 1] += m*(1 - fraction);
		QIn20[whole-1] += 0.5*m*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*(1-fraction);
		QOut20[whole-1] += m*0.5*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*r_inv*(1-fraction);
		QIn21_I[whole-1] -= m*A21*diff.y*diff.z*(1 - fraction);
		QIn21_real[whole-1] -= m*A21*diff.z*diff.x*(1 - fraction);
		QIn22_real[whole-1] += m*A22*(diff.x*diff.x - diff.y*diff.y)*(1 - fraction);
		QIn22_I[whole-1] += 2 * m*A22*(diff.x*diff.y)*(1 - fraction);
		QOut21_I[whole-1] -= m*A21*diff.y*diff.z*r_inv*(1 - fraction);
		QOut21_real[whole-1] -= m*A21*diff.z*diff.x*r_inv*(1 - fraction);
		QOut22_real[whole-1] += m*A22*(diff.x*diff.x - diff.y*diff.y)*r_inv*(1 - fraction);
		QOut22_I[whole-1] += 2 * m*A22*(diff.x*diff.y)*r_inv*(1 - fraction);*/
	}
	++resolution;
	// Sum up all the mass and moments
#ifdef RICH_MPI
	vector<double> temp(m_radius.size(), 0);
	MPI_Allreduce(&m_radius[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	m_radius = temp;
	MPI_Allreduce(&QIn20[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	QIn20 = temp;
	MPI_Allreduce(&QIn21_I[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	QIn21_I = temp;
	MPI_Allreduce(&QIn21_real[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	QIn21_real = temp;
	MPI_Allreduce(&QOut20[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	QOut20 = temp;
	MPI_Allreduce(&QOut21_I[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	QOut21_I = temp;
	MPI_Allreduce(&QOut21_real[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	QOut21_real = temp;
	MPI_Allreduce(&QOut22_real[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	QOut22_real = temp;
	MPI_Allreduce(&QOut22_I[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	QOut22_I = temp;
	MPI_Allreduce(&QIn22_I[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	QIn22_I = temp;
	MPI_Allreduce(&QIn22_real[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	QIn22_real = temp;
#else
	vector<double> temp;
#endif
	for (size_t i = 1; i < resolution; ++i)
		m_radius[i] += m_radius[i - 1];
	
	for (size_t i = 1; i < resolution; ++i)
		QIn20[i] += QIn20[i - 1];

	for (size_t i = 1; i < resolution; ++i)
		QIn21_I[i] += QIn21_I[i - 1];

	for (size_t i = 1; i < resolution; ++i)
		QIn21_real[i] += QIn21_real[i - 1];

	for (size_t i = 1; i < resolution; ++i)
		QIn22_I[i] += QIn22_I[i - 1];

	for (size_t i = 1; i < resolution; ++i)
		QIn22_real[i] += QIn22_real[i - 1];

	for (size_t i = 1; i < resolution; ++i)
		QOut20[resolution -i -1] += QOut20[resolution-i];

	for (size_t i = 1; i < resolution; ++i)
		QOut21_I[resolution - i - 1] += QOut21_I[resolution - i];

	for (size_t i = 1; i < resolution; ++i)
		QOut21_real[resolution - i - 1] += QOut21_real[resolution - i];

	for (size_t i = 1; i < resolution; ++i)
		QOut22_I[resolution - i - 1] += QOut22_I[resolution - i];

	for (size_t i = 1; i < resolution; ++i)
		QOut22_real[resolution - i - 1] += QOut22_real[resolution - i];

	// Calc acceleration
	acc.resize(N);

	for (size_t i = 0; i < N; ++i)
	{
		Vector3D point = tess.GetMeshPoint(i) - CMtot;
		double r = abs(point);
		double m = LinearInterpolation(edges_, m_radius, r);
		if (r > (2*smoothlength_))
		{
			acc[i] = point*m*(-1.0 / (r*r*r));
			Vector3D qadd = CalcQuadGrad(edges_, QIn20, QIn21_real, QIn22_real, QOut20, QOut21_real, QOut22_real, QIn21_I,
				QIn22_I, QOut21_I, QOut22_I, point);
			acc[i] -= qadd;
		}
		else
		{
			acc[i] = point*m*(-1.0 / (r*r*r + smoothlength_*smoothlength_*smoothlength_));
			Vector3D qadd = CalcQuadGrad(edges_, QIn20, QIn21_real, QIn22_real, QOut20, QOut21_real, QOut22_real, QIn21_I,
				QIn22_I, QOut21_I, QOut22_I, point);
			acc[i] -= qadd;
		}
	}

	if (output_)
	{
		size_t Np = edges_.size();
		potential.resize(N);
		vector<double> phis(Np, 0);
		phis.back() = -Mtot / edges_.back();
		for (size_t i = 1; i < Np-1; ++i)
			phis[Np - i - 1] = phis[Np - i] + m_radius[Np - i] * (1.0/edges_[Np - i] - 1.0/edges_[Np-i-1]);
		phis[0] = phis[1];
		for (size_t i = 0; i < N; ++i)
		{
			Vector3D point = tess.GetMeshPoint(i) - CMtot;
			double r = abs(point);
			size_t index = 0;
			double phi_m = LinearInterp(edges_, phis, r,index);
			double Phi = -LinearInterp2(edges_, QIn20, r,index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y) / (r*r*r*r*r);
			Phi -= LinearInterp2(edges_, QOut20, r, index)*(2 * point.z*point.z - point.x * point.x - point.y * point.y);
			Phi += 2 * LinearInterp2(edges_,QIn21_real, r, index)*point.x*point.z / (r*r*r*r*r);
			Phi += 2 * LinearInterp2(edges_,QIn21_I, r, index)*point.y*point.z / (r*r*r*r*r);
			Phi += 2 * LinearInterp2(edges_,QOut21_real, r, index)*point.x*point.z;
			Phi += 2 * LinearInterp2(edges_,QOut21_I, r, index)*point.y*point.z;
			Phi -= 2 * LinearInterp2(edges_,QIn22_real, r, index)*(point.x*point.x - point.y*point.y) / (r*r*r*r*r);
			Phi -= 4 * LinearInterp2(edges_,QIn22_I , r, index)*point.y*point.x / (r*r*r*r*r);
			Phi -= 2 * LinearInterp2(edges_,QOut22_real, r, index)*(point.x*point.x - point.y*point.y);
			Phi -= 4 * LinearInterp2(edges_,QOut22_I, r, index)*point.y*point.x;
			potential[i] = phi_m + Phi;
		}
	}
}