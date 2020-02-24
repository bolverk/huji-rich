#include "SetLoad3D.hpp"
#include <iostream>

#ifdef RICH_MPI
void SetLoad(Voronoi3D &tproc, vector<Vector3D> &points,size_t Niter, double speed, int mode, double round, bool display)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ConstNumberPerProc3D procmove(speed, round, mode);
	Voronoi3D local(tproc.GetBoxCoordinates().first, tproc.GetBoxCoordinates().second);
	local.Norg_ = points.size();
	local.del_.points_ = points;
	vector<size_t> selfindex;
	vector<vector<size_t> > sentpoints;
	vector<int> sentproc;
	int ntotal;
	for (size_t i = 0; i < Niter; ++i)
	{	
		if (display)
		{
			double load = procmove.GetLoadImbalance(local, ntotal);
			if (rank == 0)
				std::cout << "iteration " << i << " load " << load << std::endl;;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		procmove.Update(tproc, local);
		points = local.UpdateMPIPoints(tproc, rank, local.del_.points_, selfindex, sentproc, sentpoints);
		local.Norg_ = points.size();
		local.del_.points_ = points;
	}
	int total = 0;
	double load = procmove.GetLoadImbalance(local, total);
	if (rank == 0)
		std::cout << "Done SetLoad, Load is " << load << std::endl;
}

void SetLoad(Voronoi3D &tproc, vector<Vector3D> &points,vector<ComputationalCell3D> &cells, size_t Niter, double speed, 
	int mode, double round, bool display)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ConstNumberPerProc3D procmove(speed, round, mode);
	Voronoi3D local(tproc.GetBoxCoordinates().first, tproc.GetBoxCoordinates().second);
	local.Norg_ = points.size();
	local.del_.points_ = points;
	vector<size_t> selfindex;
	vector<vector<size_t> > sentpoints;
	vector<int> sentproc;
	ComputationalCell3D example;
	int total = 0;
	for (size_t i = 0; i < Niter; ++i)
	{
		if (display)
		{
			double load = procmove.GetLoadImbalance(local, total);
			if (rank == 0)
				std::cout << "iteration " << i << " load " << load << std::endl;;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		procmove.Update(tproc, local);
		points = local.UpdateMPIPoints(tproc, rank, local.del_.points_, selfindex, sentproc, sentpoints);
		local.sentpoints_ = sentpoints;
		local.self_index_ = selfindex;
		local.sentprocs_ = sentproc;
		local.Norg_ = points.size();
		MPI_exchange_data(local, cells, false,&example);
		local.del_.points_ = points;
	}
	double load = procmove.GetLoadImbalance(local, total);
	if (rank == 0)
		std::cout << "Done SetLoad, Load is " << load << std::endl;
}

#endif