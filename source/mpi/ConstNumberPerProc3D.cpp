#include "ConstNumberPerProc3D.hpp"
#ifdef RICH_MPI
#include "../misc/serializable.hpp"
#include "mpi_commands.hpp"
#include <mpi.h>
#include "HilbertProcPositions.hpp"
#endif

#ifdef RICH_MPI
namespace
{
	Vector3D GetProcCM(Tessellation3D const& tess)
	{
		size_t Ncor = tess.GetPointNo();
		Vector3D res;
		for (size_t i = 0; i < Ncor; ++i)
			res += tess.GetMeshPoint(i);
		res *= 1.0 / static_cast<double>(std::max(static_cast<size_t>(1), Ncor));
		return res;
	}
}
#endif

ConstNumberPerProc3D::~ConstNumberPerProc3D(void) {}


ConstNumberPerProc3D::ConstNumberPerProc3D(double speed, double RoundSpeed, int mode, bool Hilbert) :
	speed_(speed), RoundSpeed_(RoundSpeed), mode_(mode), Hilbert_(Hilbert), run_counter_(100) {}

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
	// Find out how many points each proc has
	vector<int> NPerProc(static_cast<size_t>(nproc));
	int mypointnumber = static_cast<int>(tlocal.GetPointNo() + 1);

	MPI_Gather(&mypointnumber, 1, MPI_INT, &NPerProc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&NPerProc[0], nproc, MPI_INT, 0, MPI_COMM_WORLD);

	int IdealPerProc = 0;
	for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		IdealPerProc += NPerProc[i];
	IdealPerProc /= nproc;
	double load = 0;
	for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		load = std::max(load, static_cast<double>(NPerProc[i]) / static_cast<double>(IdealPerProc));

	double NewSpeed = speed_;
	if (load > 3)
		NewSpeed *= std::min(std::pow(load - 2, 3), 0.15 / speed_);

	const double d = fastabs(CM - tproc.GetMeshPoint(rank));
	double dxround = 0, dyround = 0, dzround = 0;
	if (d > 0.1*R[static_cast<size_t>(rank)])
	{
		dxround = RoundSpeed_*NewSpeed*(CM.x - point.x);
		dyround = RoundSpeed_*NewSpeed*(CM.y - point.y);
		dzround = RoundSpeed_*NewSpeed*(CM.z - point.z);
	}

	if (Hilbert_ && load > 3.75 && (run_counter_ % 100 == 0))
	{
		vector<Vector3D> res = HilbertProcPositions(tlocal);
		tproc.Build(res);
		return;
	}
	++run_counter_;
	Vector3D RankCM = GetProcCM(tlocal);
	//Vector3D RankCM = tproc.GetCellCM(static_cast<size_t>(rank));
	vector<double> tosend = list_serialize(vector<Vector3D>(1, RankCM));
	vector<double> torecv(static_cast<size_t>(nproc) * 3);
	MPI_Gather(&tosend[0], 3, MPI_DOUBLE, &torecv[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&torecv[0], nproc * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	vector<Vector3D> RankCMs = list_unserialize(torecv, RankCM);
	double MyR = R[static_cast<size_t>(rank)];
	// Move point according to density
	if (mode_ == 1 || mode_ == 3)
	{
		for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		{
			if (i == static_cast<size_t>(rank))
				continue;
			Vector3D otherpoint = RankCMs[i];
			double dist = std::sqrt((point.x - otherpoint.x)*(point.x - otherpoint.x) +
				(point.y - otherpoint.y)*(point.y - otherpoint.y) + (point.z - otherpoint.z)*(point.z - otherpoint.z)
				+ 0.5*MyR * R[i]);
			double min_volume = 0.35*std::min(tproc.GetVolume(i), 0.5*dist*dist*dist);
			double merit = static_cast<double>(NPerProc[i] - IdealPerProc) / IdealPerProc;
			merit = std::max(merit, -0.5);
			double temp = merit*(point.x - otherpoint.x) * min_volume / (dist*dist*dist);
			dx -= NewSpeed*temp;
			temp = merit*(point.y - otherpoint.y) * min_volume / (dist*dist*dist);
			dy -= NewSpeed*temp;
			temp = merit*(point.z - otherpoint.z) * min_volume / (dist*dist*dist);
			dz -= NewSpeed*temp;
		}
	}
	point = tproc.GetMeshPoint(static_cast<size_t>(rank));
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

		const double neigheps = 0.05;
		for (size_t i = 0; i < Nneigh; ++i)
		{
			if (static_cast<int>(neigh[i]) >= nproc)
				continue;
			const double dist = fastabs(point - tproc.GetMeshPoint(neigh[i]));
			const double mind = neigheps*std::min(MyR, R[neigh[i]]);
			if (dist < mind)
			{
				double speed2 = NewSpeed * 10 * (mind - dist) *MyR / (dist*dist);
				dx += speed2 * (point.x - tproc.GetMeshPoint(neigh[i]).x);
				dy += speed2 * (point.y - tproc.GetMeshPoint(neigh[i]).y);
				dz += speed2 * (point.z - tproc.GetMeshPoint(neigh[i]).z);
			}
			else
			{
				Vector3D otherpoint = RankCMs[neigh[i]];
				double merit = static_cast<double>(NPerProc[static_cast<size_t>(rank)] - NPerProc[neigh[i]]) / IdealPerProc;
				double dr = fastabs(otherpoint - point);
				dx -= NewSpeed*merit*(otherpoint.x - point.x)*MyR / dr;
				dy -= NewSpeed*merit*(otherpoint.y - point.y)*MyR / dr;
				dz -= NewSpeed*merit*(otherpoint.z - point.z)*MyR / dr;
			}
		}
	}

	old_dx += dx;
	old_dy += dy;
	old_dz += dz;

	// Get min distance to other points
	vector<size_t> neigh = tproc.GetNeighbors(rank);
	size_t Nneigh = neigh.size();
	double mind_1 = 0;
	for (size_t i = 0; i < Nneigh; ++i)
	{
		if (static_cast<int>(neigh[i]) >= nproc)
			continue;
		const double dist = fastabs(point - tproc.GetMeshPoint(neigh[i]));
		mind_1 = std::max(mind_1, 1.0 / dist);
	}
	double maxR = std::min(MyR, 1.5 / mind_1);
	double r_dx = fastabs(Vector3D(old_dx, old_dy, old_dz));
	if (r_dx > NewSpeed*maxR)
	{
		old_dx *= NewSpeed*maxR / r_dx;
		old_dy *= NewSpeed*maxR / r_dx;
		old_dz *= NewSpeed*maxR / r_dx;
	}
	// Add the round cells
	double round_reduce = std::min(1.0, maxR / MyR);
	old_dx += dxround*round_reduce;
	old_dy += dyround*round_reduce;
	old_dz += dzround*round_reduce;
	r_dx = fastabs(Vector3D(old_dx, old_dy, old_dz));
	if (r_dx > NewSpeed*maxR)
	{
		old_dx *= NewSpeed*maxR / r_dx;
		old_dy *= NewSpeed*maxR / r_dx;
		old_dz *= NewSpeed*maxR / r_dx;
	}
	dx = old_dx;
	dy = old_dy;
	dz = old_dz;
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
	Vector3D cor = tproc.GetMeshPoint(static_cast<size_t>(rank)) + Vector3D(dx, dy, dz);
	// Have all processors have the same points
	tosend = list_serialize(vector<Vector3D>(1, cor));
	MPI_Gather(&tosend[0], 3, MPI_DOUBLE, &torecv[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&torecv[0], nproc * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	vector<Vector3D> cortemp = list_unserialize(torecv, cor);
	tproc.Build(cortemp);
}
#endif
