#include "SetLoad.hpp"

#ifdef RICH_MPI

void SetLoad(Tessellation &tproc,vector<Vector2D> const& points,OuterBoundary const&
	outer,int Nbest,int Niter)
{
	int rank,ws;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&ws);
	ConstNumberPerProc procmove(outer,Nbest,0.75);
	VoronoiMesh local(points,tproc,outer);
	for(int i=0;i<Niter;++i)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		procmove.Update(tproc,local);
		vector<Vector2D> cp=local.GetMeshPoints();
		cp.resize(local.GetPointNo());
		MPI_Barrier(MPI_COMM_WORLD);
		try
		{
			local.Update(cp,tproc);
		}
		catch(UniversalError const& eo)
		{
			DisplayError(eo);
		}
	}
}

#endif