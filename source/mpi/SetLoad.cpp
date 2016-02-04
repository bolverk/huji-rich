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
		outer,int Niter,double tload,double speed,int mode)
	{
		double BestLoad=100;
		vector<Vector2D> BestProc,BestMesh;
		int ws;
		MPI_Comm_size(MPI_COMM_WORLD, &ws);
		double round;
		if(mode==1)
			round=2;
		else
			round=1.6;
		ConstNumberPerProc procmove(outer,speed,round,mode);
		VoronoiMesh local(tproc, points, outer);
		vector<double> loads;
		for(int i=0;i<Niter;++i)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			//			if(rank==0)
			//				cout<<"h1"<<endl;
			procmove.Update(tproc,local);
			//			voronoi_loggers::BinLogger log("vtemp.bin");
			//			log.output(tproc);
			vector<Vector2D> cp=local.GetMeshPoints();
			cp.resize(static_cast<size_t>(local.GetPointNo()));

			/*			if(rank==0)
			{
			BestProc=tproc.GetMeshPoints();
			BestProc.resize(ws);
			WriteVector2DToFile(BestProc,"procmesh.bin");
			cout<<"h2"<<endl;
			}
			WriteVector2DToFile(cp,"localmesh"+int2str(rank)+".bin");
			*/
			try
			{
				MPI_Barrier(MPI_COMM_WORLD);
				local.Update(cp,tproc);
				//				if(rank==0)
				//				cout<<"h3"<<endl;
			}
			catch(UniversalError const& eo)
			{
				DisplayError(eo);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			double load=GetLoad(local);
			MPI_Bcast(&load, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			loads.push_back(load);
			if(load<BestLoad)
			{
				BestLoad=load;
				BestMesh=cp;
				BestProc=tproc.GetMeshPoints();
				BestProc.resize(static_cast<size_t>(ws));
			}
			if(load<tload)
			{
				break;
			}
		}
		if(mode==1)
			round=2.1;
		else
			round=1.75;
		ConstNumberPerProc procmove2(outer,speed,round,mode);
		tproc.Update(BestProc);
		MPI_Barrier(MPI_COMM_WORLD);
		local.Update(BestMesh,tproc);
		for(int i=0;i<3;++i)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			procmove2.Update(tproc,local);
			vector<Vector2D> cp=local.GetMeshPoints();
			cp.resize(static_cast<size_t>(local.GetPointNo()));
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
		points.resize(static_cast<size_t>(local.GetPointNo()));
	}
}

void SetLoad(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
	outer,int Niter,double speed,int mode)
{
	RunLoadBalance(tproc,points,outer,Niter,0,speed,mode);
}

void SetLoad(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
	outer,double TargetLoad,double speed,int mode)
{
	RunLoadBalance(tproc,points,outer,1000,TargetLoad,speed,mode);
}

#endif
