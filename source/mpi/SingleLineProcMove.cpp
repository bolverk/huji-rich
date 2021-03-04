#include "SingleLineProcMove.hpp"
#ifdef RICH_MPI
#include "../misc/serializable.hpp"
#include "../misc/simple_io.hpp"
#include "mpi_commands.hpp"
#include <mpi.h>
#include <random>
#endif

void SingleLineProcMove::Update
(
#ifdef RICH_MPI
 Tessellation3D& tproc,
 Tessellation3D const& tlocal
#else
 Tessellation3D& /*tproc*/,
 Tessellation3D const& /*tlocal*/
#endif
 
 ) const
{
#ifdef RICH_MPI
  int nproc = static_cast<int>(tproc.GetPointNo());
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int ntotal = 0;
  int nlocal = static_cast<int>(tlocal.GetPointNo());
  double load = GetLoadImbalance(tlocal, ntotal);
  if (load > 1.45)
    {
      std::mt19937_64 randgen;
      std::uniform_real_distribution<double> dist;
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
	  for (size_t i = 1; i < static_cast<size_t>(nproc); ++i)
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
	  newx[1] = 2 * allx[NperProc2] - newx[0];
	  for (int i = 2; i < nproc; ++i)
	    {
	      double x_edge = 0.5 * (newx[i - 1] + newx[i - 2]);
	      int points_done = static_cast<int>(std::upper_bound(allx.begin(), allx.end(), x_edge) - allx.begin());
	      NperProc2 = (ntotal - points_done) / (nproc - i);
	      double new_edge = allx[std::min(points_done + NperProc2, ntotal - 1)];
	      double new_loc = std::min(std::min(2 * new_edge - newx[i - 1], allx.back()), allx[points_done +
												static_cast<int>(1.5 * NperProc2)]);
	      newx[i] = std::min(new_loc, allx.back() * 0.99999);
	    }
	  write_vector(newx, "procx.txt", 8);
	}
      MPI_Bcast(&newx[0], nproc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      std::vector<Vector3D> cortemp = tproc.GetMeshPoints();
      cortemp.resize(nproc);
      for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
	cortemp[i] = Vector3D(newx[i], 0, 0);
      tproc.Build(cortemp);
    }
  else
    if (rank == 0)
      std::cout << "Load is " << load << " skipping proc move" << std::endl;
#endif
}
