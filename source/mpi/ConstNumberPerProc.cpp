#include "ConstNumberPerProc.hpp"
#ifdef RICH_MPI
#include "../misc/serializable.hpp"
#include <mpi.h>
#endif

ConstNumberPerProc::~ConstNumberPerProc(void) {}


ConstNumberPerProc::ConstNumberPerProc(OuterBoundary const& outer,double speed, double RoundSpeed, int mode) :
	outer_(outer), speed_(speed), RoundSpeed_(RoundSpeed),
	mode_(mode) {}

#ifdef RICH_MPI
void ConstNumberPerProc::Update(Tessellation& tproc, Tessellation const& tlocal)const
{
	int nproc = tproc.GetPointNo();
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<double> R(static_cast<size_t>(nproc));
	double dx = 0;
	double dy = 0;
	for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		R[i] = tproc.GetWidth(static_cast<int>(i));
	// Make cell rounder
	const Vector2D& CM = tproc.GetCellCM(rank);
	Vector2D point = tproc.GetMeshPoint(rank);
	const double d = abs(CM - tproc.GetMeshPoint(rank));
	double dxround = 0, dyround = 0;
	if (d > 0.1*R[static_cast<size_t>(rank)])
	{
		dxround = RoundSpeed_*speed_*(CM.x - point.x);
		dyround = RoundSpeed_*speed_*(CM.y - point.y);
	}
	point = CM;
	// Find out how many points each proc has
	vector<int> NPerProc(static_cast<size_t>(nproc));
	int mypointnumber = tlocal.GetPointNo() + 1;
	MPI_Gather(&mypointnumber, 1, MPI_INT, &NPerProc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&NPerProc[0], nproc, MPI_INT, 0, MPI_COMM_WORLD);
	int IdealPerProc = 0;
	for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		IdealPerProc += NPerProc[i];
	IdealPerProc /= nproc;
	// Move point according to density
	if (mode_ == 1 || mode_ == 3)
	{
		for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		{
			if (i == static_cast<size_t>(rank))
				continue;
			Vector2D otherpoint = tproc.GetCellCM(static_cast<int>(i));
			double dist = sqrt((point.x - otherpoint.x)*(point.x - otherpoint.x) +
				(point.y - otherpoint.y)*(point.y - otherpoint.y) + 0.5*R[static_cast<size_t>(rank)] * R[i]);
			double temp = (NPerProc[i] -IdealPerProc)*R[static_cast<size_t>(rank)] * R[i]*
				(point.x - otherpoint.x) / (pow(dist, 3)*IdealPerProc);
			dx -= temp;
			temp = (NPerProc[i] - IdealPerProc)*R[static_cast<size_t>(rank)] * R[i]*
				(point.y - otherpoint.y) / (pow(dist, 3)*IdealPerProc);
			dy -= temp;
		}
	}
	double old_dx = dx;
	double old_dy = dy;
	dx = 0;
	dy = 0;
	// Moving according to pressure
	if (mode_ == 1 || mode_ == 2)
	{
		const double neigheps = 0.2;
		vector<int> neigh = tproc.GetNeighbors(rank);
		for (size_t i = 0; i < neigh.size(); ++i)
		{
			if (neigh[i]>=nproc)
				continue;
			Vector2D otherpoint = tproc.GetMeshPoint(neigh[i]);
			//Vector2D otherpoint=tproc.GetCellCM(neigh[i]);
			point = tproc.GetMeshPoint(rank);
			const double dist = tproc.GetMeshPoint(rank).distance(tproc.GetMeshPoint(neigh[i]));
			if (dist < neigheps*std::min(R[static_cast<size_t>(rank)], R[static_cast<size_t>(neigh[i])]))
			{
				dx += neigheps*(point.x - tproc.GetMeshPoint(neigh[i]).x)*std::min(R[static_cast<size_t>(rank)], R[static_cast<size_t>(neigh[i])]) / dist;
				dy += neigheps*(point.y - tproc.GetMeshPoint(neigh[i]).y)*std::min(R[static_cast<size_t>(rank)], R[static_cast<size_t>(neigh[i])]) / dist;
			}
			else
			{
				dx -= (NPerProc[static_cast<size_t>(rank)] - NPerProc[static_cast<size_t>(neigh[i])])*(otherpoint.x - point.x)*R[static_cast<size_t>(rank)] / (IdealPerProc*dist);
				dy -= (NPerProc[static_cast<size_t>(rank)] - NPerProc[static_cast<size_t>(neigh[i])])*(otherpoint.y - point.y)*R[static_cast<size_t>(rank)] / (IdealPerProc*dist);
			}
		}
	}
	const double FarFraction = 0.65;
	old_dx = (old_dx>0) ? std::min(old_dx, FarFraction*speed_*R[static_cast<size_t>(rank)]) : -std::min(-old_dx, FarFraction*speed_*R[static_cast<size_t>(rank)]);
	old_dy = (old_dy > 0) ? std::min(old_dy, FarFraction*speed_*R[static_cast<size_t>(rank)]) : -std::min(-old_dy, FarFraction*speed_*R[static_cast<size_t>(rank)]);
	dx = (dx > 0) ? std::min(dx, speed_*R[static_cast<size_t>(rank)]) : -std::min(-dx, speed_*R[static_cast<size_t>(rank)]);
	dy = (dy > 0) ? std::min(dy, speed_*R[static_cast<size_t>(rank)]) : -std::min(-dy, speed_*R[static_cast<size_t>(rank)]);
	// Combine the two additions
	dx += old_dx;
	dy += old_dy;
	// Add the round cells
	dx += dxround;
	dy += dyround;
	// Make sure not out of bounds
	point = tproc.GetMeshPoint(rank);
	//	const double close=0.99999;
	const double close = 0.999;
	const double wx = tproc.GetWidth(rank);
	const double wy = wx;
	const Vector2D center(0.5*(outer_.GetGridBoundary(Left) + outer_.GetGridBoundary(Right)),
		0.5*(outer_.GetGridBoundary(Up) + outer_.GetGridBoundary(Down)));
	if (point.x + dx > outer_.GetGridBoundary(Right) - (1 - close)*wx) {
		if (outer_.GetGridBoundary(Right) - point.x < (1 - close)*wx)
			dx = -wx*(1 - close);
		else
			dx = 0.5*(outer_.GetGridBoundary(Right) - point.x);
	}
	if (point.x + dx < outer_.GetGridBoundary(Left) + (1 - close)*wx) {
		if (-outer_.GetGridBoundary(Left) + point.x < (1 - close)*wx)
			dx = wx*(1 - close);
		else
			dx = -0.5*(-outer_.GetGridBoundary(Left) + point.x);
	}
	if (point.y + dy > outer_.GetGridBoundary(Up) - (1 - close)*wy) {
		if (outer_.GetGridBoundary(Up) - point.y < (1 - close)*wy)
			dy = -wy*(1 - close);
		else
			dy = 0.5*(outer_.GetGridBoundary(Up) - point.y);
	}
	if (point.y + dy < outer_.GetGridBoundary(Down) + (1 - close)*wy) {
		if (-outer_.GetGridBoundary(Down) + point.y < (1 - close)*wy)
			dy = wy*(1 - close);
		else
			dy = 0.5*(outer_.GetGridBoundary(Down) - point.y);
	}
	Vector2D cor = tproc.GetMeshPoint(rank) + Vector2D(dx, dy);
	// Have all processors have the same points
	vector<double> tosend = list_serialize(vector<Vector2D>(1, cor));
	vector<double> torecv(static_cast<size_t>(nproc) * 2);
	MPI_Gather(&tosend[0], 2, MPI_DOUBLE, &torecv[0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&torecv[0], nproc * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	vector<Vector2D> cortemp = list_unserialize(torecv, cor);
	tproc.Update(cortemp);
}
#endif
