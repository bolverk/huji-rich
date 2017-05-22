#include "QuadrupoleGravity3D.hpp"
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
		if (*it == xi)
		{
			index = x.size() - 1;
			return y[static_cast<std::size_t>(it - x.begin())];
		}

		index = static_cast<size_t>(it - x.begin());
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

QuadrupoleGravity3D::QuadrupoleGravity3D(size_t res, double smoothlength) :r_mid_(vector<double>(res,-1)),
smoothlength_(smoothlength) 
{}

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
#endif
	CMtot *= (1.0 / Mtot);


	// GetCM, assume center for now
	size_t Nmid = r_mid_.size();
	vector<double> mr2(Nmid, 0), edges(Nmid + 1), Mrtot(Nmid + 2, 0),dr_all(Nmid, 0),Rall(Nmid+2,0);
	edges[0] = 0;
	edges.back() = abs(tess.GetBoxCoordinates().first - tess.GetBoxCoordinates().second);
	if (r_mid_[0] < 0)
	{
		for (size_t i = 0; i < Nmid; ++i)
			r_mid_[i] = (static_cast<double>(i) + 0.5)*edges.back() / static_cast<double>(Nmid);
	}
	for (size_t j = 0; j < 3; ++j)
	{
		for (size_t i = 1; i < Nmid; ++i)
			edges[i] = 0.5*(r_mid_[i] + r_mid_[i - 1]);
		mr2.assign(Nmid, 0);
		for (size_t i = 0; i < N; ++i)
		{
			size_t index = static_cast<size_t>(std::lower_bound(edges.begin(), edges.end(), abs(tess.GetMeshPoint(i) - CM)) - edges.begin() - 1);
			mr2[index] += tess.GetVolume(i)*cells[i].density;
		}
#ifdef RICH_MPI
		vector<double> temp(Nmid, 0);
		MPI_Allreduce(&mr2[0], &temp[0], static_cast<int>(Nmid), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		mr2 = temp;
#endif
		for (size_t i = 0; i < Nmid; ++i)
			mr2[i] /= (r_mid_[i] * r_mid_[i]+ smoothlength_*smoothlength_);
		for (size_t i = 1; i < Nmid + 1; ++i)
		{
			Mrtot[i] = Mrtot[i - 1] + mr2[i - 1];
			if (Mrtot[i] < (1 + 1e-7)*Mrtot[i - 1])
				Mrtot[i] = (1 + 1e-7)*Mrtot[i - 1];
		}
		Mrtot.back() = (1 + 1e-7)*Mrtot[Nmid];
		double dr = Mrtot[Nmid]/static_cast<double>(Nmid);
		Rall[0] = 0;
		for (size_t i = 0; i < Nmid; ++i)
		{
			dr_all[i] = static_cast<double>(i + 1)*dr;
			Rall[i + 1] = r_mid_[i];
		}
		Rall.back() = edges.back();
		for (size_t i = 0; i < Nmid; ++i)
			r_mid_[i] = LinearInterpolation(Mrtot, Rall, dr_all[i]);
	}
	size_t resolution = r_mid_.size();
	vector<double> m_radius(resolution), QIn20(resolution), QIn21_real(resolution),
		QIn21_I(resolution), QOut20(resolution), QOut21_real(resolution), QOut21_I(resolution),
		QOut22_I(resolution), QOut22_real(resolution), QIn22_I(resolution), QIn22_real(resolution);
	double A21 = 3.0 / 2.0, A22 = 0.25*A21;
	// Find mass and Q moments in radius
	Vector3D diff;
	for (size_t i = 0; i < N; ++i)
	{
		diff = tess.GetMeshPoint(i);
		diff -= CMtot;
		double m = tess.GetVolume(i)*cells[i].density;
		double r = abs(diff);
		if (r < r_mid_[0])
		{
			double r_inv = 1.0 / (r*r);
			double rmid = r_mid_[0];
			double r2 = rmid*rmid;
			double r_inv2= 1.0 / (rmid*rmid*rmid);
			m_radius[0] += m;
			QIn20[0] += 0.25*m*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*r_inv*r2;
			QOut20[0] += m*0.25*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*r_inv2*r_inv;
			QIn21_I[0] -= m*A21*diff.y*diff.z*r_inv*r2;
			QIn21_real[0] -= m*A21*diff.z*diff.x*r_inv*r2;
			QIn22_real[0] += m*A22*(diff.x*diff.x - diff.y*diff.y)*r_inv*r2;
			QIn22_I[0] += 2 * m*A22*(diff.x*diff.y)*r_inv*r2;
			QOut21_I[0] -= m*A21*diff.y*diff.z*r_inv2*r_inv;
			QOut21_real[0] -= m*A21*diff.z*diff.x*r_inv2*r_inv;
			QOut22_real[0] += m*A22*(diff.x*diff.x - diff.y*diff.y)*r_inv2*r_inv;
			QOut22_I[0] += 2 * m*A22*(diff.x*diff.y)*r_inv2*r_inv;
		}
		else
		{
			if ((r_mid_.back() - r) < 0)
			{
				double r_inv = 1.0 / (r*r);
				double rmid = r_mid_.back();
				double r2 = rmid*rmid;
				double r_inv2 = 1.0 / (rmid*rmid*rmid);
				m_radius[resolution-1] += m;
				QIn20[resolution - 1] += 0.25*m*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*r_inv*r2;
				QOut20[resolution - 1] += m*0.25*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*r_inv2*r_inv;
				QIn21_I[resolution - 1] -= m*A21*diff.y*diff.z*r_inv*r2;
				QIn21_real[resolution - 1] -= m*A21*diff.z*diff.x*r_inv*r2;
				QIn22_real[resolution - 1] += m*A22*(diff.x*diff.x - diff.y*diff.y)*r_inv*r2;
				QIn22_I[resolution - 1] += 2 * m*A22*(diff.x*diff.y)*r_inv*r2;
				QOut21_I[resolution - 1] -= m*A21*diff.y*diff.z*r_inv2*r_inv;
				QOut21_real[resolution - 1] -= m*A21*diff.z*diff.x*r_inv2*r_inv;
				QOut22_real[resolution - 1] += m*A22*(diff.x*diff.x - diff.y*diff.y)*r_inv2*r_inv;
				QOut22_I[resolution - 1] += 2 * m*A22*(diff.x*diff.y)*r_inv2*r_inv;
			}
			else
			{
				size_t whole = static_cast<size_t>(std::lower_bound(r_mid_.begin(), r_mid_.end(), r) - r_mid_.begin());
				if (whole > 0)
					--whole;
				double rmid0 = r_mid_[whole];
				double rmid1 = r_mid_[whole + 1];
				double fraction = 1 - (r - rmid0) / (rmid1-rmid0);
				double r_inv = 1.0 / (r*r);
				double r_invmid0 = 1.0 / (rmid0*rmid0*rmid0);
				double r20 = rmid0*rmid0;
				double r_invmid1 = 1.0 / (rmid1*rmid1*rmid1);
				double r21 = rmid1*rmid1;
				m_radius[whole] += m*fraction;
				QIn20[whole] += 0.25*m*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*r_inv*r20*fraction;
				QOut20[whole] += m*0.25*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*r_inv*r_invmid0*fraction;
				QIn21_I[whole] -= m*A21*diff.y*diff.z*r_inv*r20*fraction;
				QIn21_real[whole] -= m*A21*diff.z*diff.x*r_inv*r20*fraction;
				QIn22_real[whole] += m*A22*(diff.x*diff.x - diff.y*diff.y)*r_inv*r20*fraction;
				QIn22_I[whole] += 2 * m*A22*(diff.x*diff.y)*r_inv*r20*fraction;
				QOut21_I[whole] -= m*A21*diff.y*diff.z*r_inv*r_invmid0*fraction;
				QOut21_real[whole] -= m*A21*diff.z*diff.x*r_inv*r_invmid0*fraction;
				QOut22_real[whole] += m*A22*(diff.x*diff.x - diff.y*diff.y)*r_inv*r_invmid0*fraction;
				QOut22_I[whole] += 2 * m*A22*(diff.x*diff.y)*r_inv*r_invmid0*fraction;

				m_radius[whole + 1] += m*(1 - fraction);
				QIn20[whole+1] += 0.5*m*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*r_inv*r21*(1-fraction);
				QOut20[whole+1] += m*0.5*(2 * diff.z*diff.z - diff.x*diff.x - diff.y*diff.y)*r_inv*r_invmid1*(1-fraction);
				QIn21_I[whole+1] -= m*A21*diff.y*diff.z*r_inv*r21*(1 - fraction);
				QIn21_real[whole+1] -= m*A21*diff.z*diff.x*r_inv*r21*(1 - fraction);
				QIn22_real[whole+1] += m*A22*(diff.x*diff.x - diff.y*diff.y)*r_inv*r21*(1 - fraction);
				QIn22_I[whole+1] += 2 * m*A22*(diff.x*diff.y)*r_inv*r21*(1 - fraction);
				QOut21_I[whole+1] -= m*A21*diff.y*diff.z*r_inv*r_invmid1*(1 - fraction);
				QOut21_real[whole+1] -= m*A21*diff.z*diff.x*r_inv*r_invmid1*(1 - fraction);
				QOut22_real[whole+1] += m*A22*(diff.x*diff.x - diff.y*diff.y)*r_inv*r_invmid1*(1 - fraction);
				QOut22_I[whole+1] += 2 * m*A22*(diff.x*diff.y)*r_inv*r_invmid1*(1 - fraction);
			}
		}
	}
	// Sum up all the mass and moments
#ifdef RICH_MPI
	vector<double> temp(m_radius.size(), 0);
	MPI_Allreduce(&m_radius[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	m_radius = temp;
	MPI_Allreduce(&QIn20[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	QIn20 = temp;
	MPI_Allreduce(&QIn21_I[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	QIn21_I = temp;
	MPI_Allreduce(&QIn21_real[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	QIn21_real = temp;
	MPI_Allreduce(&QOut20[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	QOut20 = temp;
	MPI_Allreduce(&QOut21_I[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	QOut21_I = temp;
	MPI_Allreduce(&QOut21_real[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	QOut21_real = temp;
	MPI_Allreduce(&QOut22_real[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	QOut22_real = temp;
	MPI_Allreduce(&QOut22_I[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	QOut22_I = temp;
	MPI_Allreduce(&QIn22_I[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	QIn22_I = temp;
	MPI_Allreduce(&QIn22_real[0], &temp[0], static_cast<int>(resolution), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	QIn22_real = temp;
#else
	vector<double> temp;
#endif
	temp = m_radius;
	for (size_t i = 1; i < resolution; ++i)
		m_radius[i] += m_radius[i - 1];
	m_radius[0] *= 0.5;
	for (size_t i = 1; i < resolution; ++i)
		m_radius[i] -= temp[i] * 0.5;
	
	temp = QIn20;
	for (size_t i = 1; i < resolution; ++i)
		QIn20[i] += QIn20[i - 1];
	QIn20[0] *= 0.5;
	for (size_t i = 1; i < resolution; ++i)
		QIn20[i] -= temp[i] * 0.5;

	temp = QIn21_I;
	for (size_t i = 1; i < resolution; ++i)
		QIn21_I[i] += QIn21_I[i - 1];
	QIn21_I[0] *= 0.5;
	for (size_t i = 1; i < resolution; ++i)
		QIn21_I[i] -= temp[i] * 0.5;

	temp = QIn21_real;
	for (size_t i = 1; i < resolution; ++i)
		QIn21_real[i] += QIn21_real[i - 1];
	QIn21_real[0] *= 0.5;
	for (size_t i = 1; i < resolution; ++i)
		QIn21_real[i] -= temp[i] * 0.5;

	temp = QIn22_I;
	for (size_t i = 1; i < resolution; ++i)
		QIn22_I[i] += QIn22_I[i - 1];
	QIn22_I[0] *= 0.5;
	for (size_t i = 1; i < resolution; ++i)
		QIn22_I[i] -= temp[i] * 0.5;

	temp = QIn22_real;
	for (size_t i = 1; i < resolution; ++i)
		QIn22_real[i] += QIn22_real[i - 1];
	QIn22_real[0] *= 0.5;
	for (size_t i = 1; i < resolution; ++i)
		QIn22_real[i] -= temp[i] * 0.5;

	temp = QOut20;
	for (size_t i = 1; i < resolution; ++i)
		QOut20[resolution -i -1] += QOut20[resolution-i];
	QOut20.back() *= 0.5;
	for (size_t i = 0; i < resolution-1; ++i)
		QOut20[resolution - i - 1] -= temp[resolution - i - 1] * 0.5;

	temp = QOut21_I;
	for (size_t i = 1; i < resolution; ++i)
		QOut21_I[resolution - i - 1] += QOut21_I[resolution - i];
	QOut21_I.back() *= 0.5;
	for (size_t i = 0; i < resolution - 1; ++i)
		QOut21_I[resolution - i - 1] -= temp[resolution - i - 1] * 0.5;

	temp = QOut21_real;
	for (size_t i = 1; i < resolution; ++i)
		QOut21_real[resolution - i - 1] += QOut21_real[resolution - i];
	QOut21_real.back() *= 0.5;
	for (size_t i = 0; i < resolution - 1; ++i)
		QOut21_real[resolution - i - 1] -= temp[resolution - i - 1] * 0.5;

	temp = QOut22_I;
	for (size_t i = 1; i < resolution; ++i)
		QOut22_I[resolution - i - 1] += QOut22_I[resolution - i];
	QOut22_I.back() *= 0.5;
	for (size_t i = 0; i < resolution - 1; ++i)
		QOut22_I[resolution - i - 1] -= temp[resolution - i - 1] * 0.5;

	temp = QOut22_real;
	for (size_t i = 1; i < resolution; ++i)
		QOut22_real[resolution - i - 1] += QOut22_real[resolution - i];
	QOut22_real.back() *= 0.5;
	for (size_t i = 0; i < resolution - 1; ++i)
		QOut22_real[resolution - i - 1] -= temp[resolution - i - 1] * 0.5;

	// Calc acceleration
	vector<double> r_list(r_mid_);
	r_list.insert(r_list.begin(), 0);
	r_list.push_back(r_list.back() + 0.5*r_list.back()-r_list[resolution-1]);
	QIn20.insert(QIn20.begin(), QIn20[0]);
	QIn20.push_back(QIn20.back());
	QOut20.insert(QOut20.begin(), QOut20[0]);
	QOut20.push_back(QOut20.back());
	QIn21_I.insert(QIn21_I.begin(), QIn21_I[0]);
	QIn21_I.push_back(QIn21_I.back());
	QIn21_real.insert(QIn21_real.begin(), QIn21_real[0]);
	QIn21_real.push_back(QIn21_real.back());
	QIn22_I.insert(QIn22_I.begin(), QIn22_I[0]);
	QIn22_I.push_back(QIn22_I.back());
	QIn22_real.insert(QIn22_real.begin(), QIn22_real[0]);
	QIn22_real.push_back(QIn22_real.back());
	QOut21_I.insert(QOut21_I.begin(), QOut21_I[0]);
	QOut21_I.push_back(QOut21_I.back());
	QOut21_real.insert(QOut21_real.begin(), QOut21_real[0]);
	QOut21_real.push_back(QOut21_real.back());
	QOut22_I.insert(QOut22_I.begin(), QOut22_I[0]);
	QOut22_I.push_back(QOut22_I.back());
	QOut22_real.insert(QOut22_real.begin(), QOut22_real[0]);
	QOut22_real.push_back(QOut22_real.back());
	acc.resize(N);
	for (size_t i = 0; i < N; ++i)
	{
		Vector3D point = tess.GetMeshPoint(i) - CMtot;
		double r = abs(point);
		double m = LinearInterpolation(r_list, m_radius, r);
		if (r > (2*smoothlength_))
		{
			acc[i] = point*m*(-1.0 / (r*r*r));
			Vector3D qadd = CalcQuadGrad(r_list, QIn20, QIn21_real, QIn22_real, QOut20, QOut21_real, QOut22_real, QIn21_I,
				QIn22_I, QOut21_I, QOut22_I, point);
			acc[i] -= qadd;
		}
		else
			acc[i] = point*m*(-1.0 / (r*r*r+smoothlength_*smoothlength_*smoothlength_));
	}
}