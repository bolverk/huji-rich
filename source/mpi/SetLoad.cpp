#include "SetLoad.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif


#ifdef RICH_MPI

namespace
{
	double GetLoad(Tessellation const& tess)
	{
		int ws;
		MPI_Comm_size(MPI_COMM_WORLD, &ws);
		int n=tess.GetPointNo();
		int result=10;
		int total=n;
		MPI_Reduce(&n, &result, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&n, &total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		return result*ws*1.0/(total*1.0);
	}

	void RunLoadBalance(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
		outer,int Niter,double tload,double speed,int mode,bool Rmin)
	{
		int ws;
		MPI_Comm_size(MPI_COMM_WORLD, &ws);
		int rank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		double round;
		if(mode==1)
			round=2;
		else
			round=1.6;
		vector<size_t> selfindex;
		vector<vector<int> > sentpoints;
		vector<int> sentproc;
		ConstNumberPerProc procmove(outer,speed,round,mode,Rmin);
		VoronoiMesh local;
		local.GetMeshPoints() = points;
		local.SetPointNo(static_cast<int>(points.size()));
		for(int i=0;i<Niter;++i)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			procmove.Update(tproc,local);
			points = local.UpdateMPIPoints(tproc, rank, points, &outer, selfindex, sentproc, sentpoints);
			local.GetMeshPoints() = points;
			local.SetPointNo(static_cast<int>(points.size()));
			MPI_Barrier(MPI_COMM_WORLD);
			double load=GetLoad(local);
			MPI_Bcast(&load, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if(load<tload)
			{
				break;
			}
		}
	}
}

void SetLoad(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
	outer,int Niter,double speed,int mode,bool Rmin)
{
	RunLoadBalance(tproc,points,outer,Niter,0,speed,mode,Rmin);
}

void SetLoad(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
	outer,double TargetLoad,double speed,int mode,bool Rmin)
{
	RunLoadBalance(tproc,points,outer,1000,TargetLoad,speed,mode,Rmin);
}

#endif
