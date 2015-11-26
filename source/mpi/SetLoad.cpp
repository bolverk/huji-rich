#include "SetLoad.hpp"
#ifdef RICH_MPI
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#endif


#ifdef RICH_MPI

namespace
{
	double GetLoad(Tessellation const& tess)
	{
		const boost::mpi::communicator world;
		const int ws = world.size();
		int n=tess.GetPointNo();
		int result=10;
		int total=n;
		boost::mpi::reduce(world, n, result,boost::mpi::maximum<int>(), 0);
		boost::mpi::reduce(world, n, total, std::plus<int>(), 0);
		return result*ws*1.0/(total*1.0);
	}

	void RunLoadBalance(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
		outer,int Nbest,int Niter,double tload,double speed,int mode)
	{
		double BestLoad=100;
		vector<Vector2D> BestProc,BestMesh;
		const boost::mpi::communicator world;
		const int ws = world.size();
		double round;
		if(mode==1)
			round=2;
		else
			round=1.6;
		ConstNumberPerProc procmove(outer,Nbest,speed,round,mode);
		VoronoiMesh local(tproc, points, outer);
		vector<double> loads;
		for(int i=0;i<Niter;++i)
		{
			world.barrier();
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
				local.Update(cp,tproc);
				//				if(rank==0)
				//				cout<<"h3"<<endl;
			}
			catch(UniversalError const& eo)
			{
				DisplayError(eo);
			}
			world.barrier();
			double load=GetLoad(local);
			boost::mpi::broadcast(world, load, 0);
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
		ConstNumberPerProc procmove2(outer,Nbest,speed,round,mode);
		tproc.Update(BestProc);
		world.barrier();
		local.Update(BestMesh,tproc);
		for(int i=0;i<3;++i)
		{
			world.barrier();
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
	outer,int Nbest,int Niter,double speed,int mode)
{
	RunLoadBalance(tproc,points,outer,Nbest,Niter,0,speed,mode);
}

void SetLoad(Tessellation &tproc,vector<Vector2D> &points,OuterBoundary const&
	outer,int Nbest,double TargetLoad,double speed,int mode)
{
	RunLoadBalance(tproc,points,outer,Nbest,1000,TargetLoad,speed,mode);
}

#endif
