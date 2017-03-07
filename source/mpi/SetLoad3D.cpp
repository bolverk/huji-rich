#include "SetLoad3D.hpp"

#ifdef RICH_MPI
void SetLoad(Voronoi3D &tproc, vector<Vector3D> &points,size_t Niter, double speed, int mode)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<Vector2D> BestProc, BestMesh;
	double round;
	if (mode == 1)
		round = 2;
	else
		round = 1.6;
	ConstNumberPerProc3D procmove(speed, round, mode);
	Voronoi3D local(tproc.GetBoxCoordinates().first, tproc.GetBoxCoordinates().second);
	local.Norg_ = points.size();
	local.del_.points_ = points;
	vector<size_t> selfindex;
	vector<vector<size_t> > sentpoints;
	vector<int> sentproc;
	for (size_t i = 0; i < Niter; ++i)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		procmove.Update(tproc, local);
		points = local.UpdateMPIPoints(tproc, rank, local.del_.points_, selfindex, sentproc, sentpoints);
		local.Norg_ = points.size();
		local.del_.points_ = points;
	}
}
#endif