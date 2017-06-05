#include "ConstNumberPerProc3D.hpp"
#ifdef RICH_MPI
#include "../misc/serializable.hpp"
#include "mpi_commands.hpp"
#include <mpi.h>
#endif

ConstNumberPerProc3D::~ConstNumberPerProc3D(void) {}


ConstNumberPerProc3D::ConstNumberPerProc3D(double speed, double RoundSpeed, int mode) :
	speed_(speed), RoundSpeed_(RoundSpeed),	mode_(mode) {}

#ifdef RICH_MPI
void ConstNumberPerProc3D::Update(Tessellation3D& tproc, Tessellation3D const& tlocal)const
{
	int nproc = static_cast<int>(tproc.GetPointNo());
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<double> R(static_cast<size_t>(nproc));
	double dx = 0;
	double dy = 0;
	double dz = 0;
	for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		R[i] = tproc.GetWidth(i);
	// Make cell rounder
	const Vector3D& CM = tproc.GetCellCM(rank);
	Vector3D point = tproc.GetMeshPoint(rank);
	const double d = abs(CM - tproc.GetMeshPoint(rank));
	double dxround = 0, dyround = 0, dzround = 0;
	if (d > 0.1*R[static_cast<size_t>(rank)])
	{
		dxround = RoundSpeed_*speed_*(CM.x - point.x);
		dyround = RoundSpeed_*speed_*(CM.y - point.y);
		dzround = RoundSpeed_*speed_*(CM.z - point.z);
	}
	point = CM;
	// Find out how many points each proc has
	vector<int> NPerProc(static_cast<size_t>(nproc));
	int mypointnumber = static_cast<int>(tlocal.GetPointNo() + 1);
	MPI_Gather(&mypointnumber, 1, MPI_INT, &NPerProc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&NPerProc[0], nproc, MPI_INT, 0, MPI_COMM_WORLD);
	int IdealPerProc = 0;
	for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		IdealPerProc += NPerProc[i];
	IdealPerProc /= nproc;
	double load = 0;
	for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		load = std::max(load, static_cast<double>(NPerProc[i]) / static_cast<double>(IdealPerProc));
	// Move point according to density
	if (mode_ == 1 || mode_ == 3)
	{
		for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		{
			if (i == static_cast<size_t>(rank))
				continue;
			Vector3D otherpoint = tproc.GetCellCM(i);
			double dist = sqrt((point.x - otherpoint.x)*(point.x - otherpoint.x) +
				(point.y - otherpoint.y)*(point.y - otherpoint.y) + (point.z - otherpoint.z)*(point.z - otherpoint.z)
				+ 0.5*R[static_cast<size_t>(rank)] * R[i]);
			double temp = (NPerProc[i] - IdealPerProc)* 
				(point.x - otherpoint.x) / (pow(dist / R[static_cast<size_t>(rank)], 3)*IdealPerProc);
			dx -= temp;
			temp = (NPerProc[i] - IdealPerProc)* 
				(point.y - otherpoint.y) / (pow(dist/ R[static_cast<size_t>(rank)], 3)*IdealPerProc);
			dy -= temp;
			temp = (NPerProc[i] - IdealPerProc)* 
				(point.z - otherpoint.z) / (pow(dist/ R[static_cast<size_t>(rank)], 3)*IdealPerProc);
			dz -= temp;
		}
	}
	double old_dx = dx;
	double old_dy = dy;
	double old_dz = dz;
	dx = 0;
	dy = 0;
	dz = 0;
	// Moving according to pressure
	if (mode_ == 1 || mode_ == 2)
	{
		vector<size_t> neigh = tproc.GetNeighbors(rank);
		size_t Nneigh = neigh.size();
		
		const double neigheps = 0.2;
		for (size_t i = 0; i < Nneigh; ++i)
		{
			if (static_cast<int>(neigh[i]) >= nproc)
				continue;
			Vector3D otherpoint = tproc.GetMeshPoint(neigh[i]);
			point = tproc.GetMeshPoint(rank);
			const double dist = abs(tproc.GetMeshPoint(static_cast<size_t>(rank))-(tproc.GetMeshPoint(neigh[i])));
			if (dist < neigheps*std::min(R[static_cast<size_t>(rank)], R[neigh[i]]))
			{
				dx += (point.x - tproc.GetMeshPoint(neigh[i]).x)*std::min(R[static_cast<size_t>(rank)], R[neigh[i]]) / (dist*neigheps);
				dy += (point.y - tproc.GetMeshPoint(neigh[i]).y)*std::min(R[static_cast<size_t>(rank)], R[neigh[i]]) / (dist*neigheps);
				dz += (point.z - tproc.GetMeshPoint(neigh[i]).z)*std::min(R[static_cast<size_t>(rank)], R[neigh[i]]) / (dist*neigheps);
			}
			else
			{
				double merit = (NPerProc[static_cast<size_t>(rank)] - NPerProc[neigh[i]]) / IdealPerProc;
				if (merit < 0)
					merit = std::max(merit, -2.0);
				else
					merit = std::min(merit, 2.0);
				dx -= merit*(otherpoint.x - point.x)*R[static_cast<size_t>(rank)];
				dy -= merit*(otherpoint.y - point.y)*R[static_cast<size_t>(rank)];
				dz -= merit*(otherpoint.z - point.z)*R[static_cast<size_t>(rank)];
			}
		}
	}
	const double FarFraction = load > 3 ? 1.1 : 0.65;
	old_dx = (old_dx>0) ? std::min(old_dx, FarFraction*speed_*R[static_cast<size_t>(rank)]) : -std::min(-old_dx, FarFraction*speed_*R[static_cast<size_t>(rank)]);
	old_dy = (old_dy > 0) ? std::min(old_dy, FarFraction*speed_*R[static_cast<size_t>(rank)]) : -std::min(-old_dy, FarFraction*speed_*R[static_cast<size_t>(rank)]);
	old_dz = (old_dz > 0) ? std::min(old_dz, FarFraction*speed_*R[static_cast<size_t>(rank)]) : -std::min(-old_dz, FarFraction*speed_*R[static_cast<size_t>(rank)]);
	dx = (dx > 0) ? std::min(dx, speed_*R[static_cast<size_t>(rank)]) : -std::min(-dx, speed_*R[static_cast<size_t>(rank)]);
	dy = (dy > 0) ? std::min(dy, speed_*R[static_cast<size_t>(rank)]) : -std::min(-dy, speed_*R[static_cast<size_t>(rank)]);
	dz = (dz > 0) ? std::min(dz, speed_*R[static_cast<size_t>(rank)]) : -std::min(-dz, speed_*R[static_cast<size_t>(rank)]);
	// Combine the two additions
	dx += old_dx;
	dy += old_dy;
	dz += old_dz;
	// Add the round cells
	dx += dxround;
	dy += dyround;
	dz += dzround;
	// Make sure not out of bounds
	point = tproc.GetMeshPoint(rank);
	const double close = 0.999;
	const double wx = tproc.GetWidth(rank);
	Vector3D ll = tproc.GetBoxCoordinates().first;
	Vector3D ur = tproc.GetBoxCoordinates().second;
	const Vector3D center(0.5*(ll + ur));
	if (point.x + dx > (ur.x - (1 - close)*wx)) 
	{
		if ((ur.x - point.x) < ((1 - close)*wx))
			dx = -wx*(1 - close);
		else
			dx = 0.5*(ur.x - point.x);
	}
	if ((point.x + dx) < (ll.x + (1 - close)*wx)) 
	{
		if ((-ll.x + point.x) < ((1 - close)*wx))
			dx = wx*(1 - close);
		else
			dx = -0.5*(-ll.x + point.x);
	}
	if (point.y + dy > (ur.y - (1 - close)*wx))
	{
		if ((ur.y - point.y) < ((1 - close)*wx))
			dy = -wx*(1 - close);
		else
			dy = 0.5*(ur.y - point.y);
	}
	if ((point.y + dy) < (ll.y + (1 - close)*wx))
	{
		if ((-ll.y + point.y) < ((1 - close)*wx))
			dy = wx*(1 - close);
		else
			dy = -0.5*(-ll.y + point.y);
	}
	if (point.z + dz > (ur.z - (1 - close)*wx))
	{
		if ((ur.z - point.z) < ((1 - close)*wx))
			dz = -wx*(1 - close);
		else
			dz = 0.5*(ur.z - point.z);
	}
	if ((point.z + dz) < (ll.z + (1 - close)*wx))
	{
		if ((-ll.z + point.z) < ((1 - close)*wx))
			dz = wx*(1 - close);
		else
			dz = -0.5*(-ll.z + point.z);
	}
	Vector3D cor = tproc.GetMeshPoint(static_cast<size_t>(rank)) + Vector3D(dx, dy,dz);
	// Have all processors have the same points
	vector<double> tosend = list_serialize(vector<Vector3D>(1, cor));
	vector<double> torecv(static_cast<size_t>(nproc) * 3);
	MPI_Gather(&tosend[0], 3, MPI_DOUBLE, &torecv[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&torecv[0], nproc * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	vector<Vector3D> cortemp = list_unserialize(torecv, cor);
	tproc.Build(cortemp);
}
#endif