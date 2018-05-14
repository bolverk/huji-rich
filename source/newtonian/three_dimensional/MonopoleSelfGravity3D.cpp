#include "MonopoleSelfGravity3D.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif

MonopoleSelfGravity3D::MonopoleSelfGravity3D(size_t resolution,double smoothlength) :resolution_(resolution),
smoothlength_(smoothlength){}

void MonopoleSelfGravity3D::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
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
	CMtot *= (1.0/Mtot);
	// Find max radius
	double Rmax = std::max(abs(tess.GetBoxCoordinates().first - CMtot), abs(tess.GetBoxCoordinates().second - CMtot));
	double dr = Rmax / (static_cast<double>(resolution_)-1);
	// Create radius list
	vector<double> m_radius(resolution_),r_list(resolution_);
	for (size_t i = 0; i < resolution_; ++i)
		r_list[i] = dr*static_cast<double>(i);
	// Find mass in radius
	for (size_t i = 0; i < N; ++i)
	{
		double r = fastabs(tess.GetMeshPoint(i)-CMtot);
		double fraction = r / dr;
		size_t whole = static_cast<size_t>(fraction);
		double m = tess.GetVolume(i)*cells[i].density;
		if (whole == 0)
			m_radius[0] += m;
		else
		{
			double fraction2 = fraction - static_cast<double>(whole);
			m_radius[whole] += m*(1 - fraction2);
			m_radius[whole + 1] += m*fraction2;
		}
	}
	// Sum up all the mass
#ifdef RICH_MPI
	vector<double> temp(m_radius.size(), 0);
	MPI_Allreduce(&m_radius[0], &temp[0], static_cast<int>(resolution_), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	m_radius = temp;
#endif
	for (size_t i = 1; i < resolution_; ++i)
	{
		m_radius[i] += m_radius[i - 1];
	}
	// Calc acceleration
	acc.resize(N);
	for (size_t i = 0; i < N; ++i)
	{
		Vector3D point = tess.GetMeshPoint(i)-CMtot;
		double r = fastabs(point);
		double m = LinearInterpolation(r_list, m_radius, r);
		if(r>(smoothlength_))
			acc[i] = point*m*(-1.0/(r*r*r));
		else
			acc[i] = point*m*(-1.0 / (smoothlength_*smoothlength_*smoothlength_));
	}
}