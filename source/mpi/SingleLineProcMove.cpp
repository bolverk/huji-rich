#include "SingleLineProcMove.hpp"
#ifdef RICH_MPI
#include "../misc/serializable.hpp"
#include "mpi_commands.hpp"
#include <mpi.h>
#endif

void SingleLineProcMove::Update(Tessellation3D& tproc, Tessellation3D const& tlocal) const
{
#ifdef RICH_MPI
	int nproc = static_cast<int>(tproc.GetPointNo());
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int ntotal = 0;
	int nlocal = tlocal.GetPointNo();
	double load = GetLoadImbalance(tlocal, ntotal);
	if (load > 1.55)
	{
		std::vector<double> allx, xlocal(nlocal);
		if (rank == 0)
		{
			allx.resize(ntotal);
		}
		for (int i = 0; i < nlocal; ++i)
			xlocal[i] = tlocal.GetMeshPoint(i).x;
		vector<int> NPerProc(static_cast<size_t>(nproc)), disp;
		MPI_Gather(&nlocal, 1, MPI_INT, &NPerProc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			disp.resize(nproc, 0);
			for (size_t i = 1; i < nproc; ++i)
				disp[i] = disp[i - 1] + NPerProc[i - 1];
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gatherv(&xlocal[0], nlocal, MPI_DOUBLE, &allx[0], &NPerProc[0], &disp[0],
			MPI_DOUBLE, 0, MPI_COMM_WORLD);
		std::vector<double> newx(nproc);
		if (rank == 0)
		{
			std::sort(allx.begin(), allx.end());
			int NperProc2 = ntotal / nproc + 1;
			newx[0] = allx[NperProc2] - 0.5 * (allx[2 * NperProc2] - allx[NperProc2]);
			for (int i = 1; i < nproc; ++i)
			{
				int index = std::min(NperProc2 * i, ntotal - 1);
				newx[i] = 2 * allx.at(index) - 0.001 * (allx.at(index) - allx.at(index - 1))
					- newx[i - 1];
			}
		}
		MPI_Bcast(&newx[0], nproc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		std::vector<Vector3D> cortemp = tproc.GetMeshPoints();
		cortemp.resize(nproc);
		for (size_t i = 0; i < nproc; ++i)
			cortemp[i] = Vector3D(newx[i], 0, 0);
		tproc.Build(cortemp);
	}
	else
		if (rank == 0)
			std::cout << "Load is " << load << " skipping proc move" << std::endl;
#endif
}
