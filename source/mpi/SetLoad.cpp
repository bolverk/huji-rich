#include "SetLoad.hpp"

#ifdef RICH_MPI

namespace
{
	double GetLoad(Tessellation const& tess)
	{
		int ws;
		MPI_Comm_size(MPI_COMM_WORLD,&ws);
		int n=tess.GetPointNo();
		int result=0;
		int total=n;
		MPI_Reduce(&n,&result,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
		MPI_Reduce(&n,&total,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		return result*ws*1.0/(total*1.0);
	}

	void RunLoadBalance(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
		outer,int Nbest,int Niter,double tload)
	{
		int rank,ws;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&ws);
		ConstNumberPerProc procmove(outer,Nbest,0.5);
		VoronoiMesh local(points,tproc,outer);
		for(int i=0;i<Niter;++i)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			procmove.Update(tproc,local);
			vector<Vector2D> cp=local.GetMeshPoints();
			cp.resize(local.GetPointNo());
			MPI_Barrier(MPI_COMM_WORLD);
			double load=GetLoad(local);
			if(rank==0&&i%10==0)
				cout<<"Setting the load balance, iteration="<<i<<" load="<<load<<endl;
			if(load<tload)
			{
				points=local.GetMeshPoints();
				points.resize(local.GetPointNo());
				return;
			}
			try
			{
				local.Update(cp,tproc);
			}
			catch(UniversalError const& eo)
			{
				DisplayError(eo);
			}
		}
		points=local.GetMeshPoints();
		points.resize(local.GetPointNo());
	}

}

void SetLoad(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
	outer,int Nbest,int Niter)
{
	RunLoadBalance(tproc,points,outer,Nbest,Niter,0);
}

void SetLoad(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
	outer,int Nbest,double TargetLoad)
{
	RunLoadBalance(tproc,points,outer,Nbest,1000,TargetLoad);
}

#endif