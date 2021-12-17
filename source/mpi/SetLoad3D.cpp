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

void SetLoad(HDSim3D& sim, size_t Niter, double speed, int mode,
	double round, bool display, bool const tess_rebuild)
{
	int total = 0;
	int rank = 0;
	Conserved3D edummy;
	ComputationalCell3D cdummy;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ConstNumberPerProc3D procmove(speed, round, mode);
	std::vector<size_t> selfindex;
	std::vector<vector<size_t> > sentpoints;
	std::vector<int> sentproc;
	std::vector<Vector3D> newpoints;
	Tessellation3D& tproc = sim.getProcTesselation();
	Tessellation3D& local = sim.getTesselation();
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
		std::vector<Vector3D> &oldpoints = local.accessMeshPoints();
		oldpoints.resize(local.GetPointNo());
		newpoints = local.UpdateMPIPoints(tproc, rank, oldpoints, selfindex, sentproc, sentpoints);
		MPI_Barrier(MPI_COMM_WORLD);
		local.GetPointNo() = newpoints.size();
		local.GetSentPoints() = sentpoints;
		local.GetSentProcs() = sentproc;
		local.GetSelfIndex() = selfindex;
		local.accessMeshPoints() = newpoints;
		// Keep relevant points
		MPI_exchange_data(local, sim.getExtensives(), false, &edummy);
		MPI_exchange_data(local, sim.getCells(), false, &cdummy);
	}
	if(tess_rebuild)
	{
		local.Build(newpoints, tproc);
		MPI_exchange_data(local, sim.getCells(), true, &cdummy);
	}
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
