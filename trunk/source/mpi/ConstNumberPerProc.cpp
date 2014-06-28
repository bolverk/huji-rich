#include "ConstNumberPerProc.hpp"

ConstNumberPerProc::~ConstNumberPerProc(void){}

ConstNumberPerProc::ConstNumberPerProc(OuterBoundary const& outer,int npercell,
	double speed,double RoundSpeed):
outer_(outer),PointsPerProc_(npercell),speed_(speed),RoundSpeed_(RoundSpeed){}

void ConstNumberPerProc::Update(Tessellation &tproc,Tessellation const& tlocal)const
{
#ifdef RICH_MPI
	int nproc=tproc.GetPointNo();
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	vector<double> R(nproc);
	double dx=0;
	double dy=0;
	for(int i=0;i<nproc;++i)
		R[i]=tproc.GetWidth(i);
	// Make cell rounder
	const Vector2D CM=tproc.GetCellCM(rank);
	const Vector2D point=tproc.GetMeshPoint(rank);
	const double d=abs(CM-tproc.GetMeshPoint(rank));
	double dxround=0,dyround=0;
	if(d>0.1*R[rank])
	{
		dxround=RoundSpeed_*speed_*(CM.x-point.x);
		dyround=RoundSpeed_*speed_*(CM.y-point.y);
	}
	// Find out how many points each proc has
	vector<int> NPerProc(nproc);
	int mypointnumber=tlocal.GetPointNo();
	MPI_Gather(&mypointnumber,1,MPI_INT,&NPerProc[0],1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&NPerProc[0],nproc,MPI_INT,0,MPI_COMM_WORLD);
	// Move point according to density
	/*
	for(int i=0;i<nproc;++i)
	{
		Vector2D otherpoint=tproc.GetMeshPoint(i);
		double dist=sqrt((point.x-otherpoint.x)*(point.x-otherpoint.x)+
			(point.y-otherpoint.y)*(point.y-otherpoint.y)+0.2*R[rank]*R[i]);
		double temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.x-otherpoint.x)/(dist*dist*dist);
		dx+=temp;
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.y-otherpoint.y)/(dist*dist*dist);
		dy+=temp;
		// Right side
		otherpoint.x=2*outer_.GetGridBoundary(Right)-otherpoint.x;
		dist=sqrt((point.x-otherpoint.x)*(point.x-otherpoint.x)+
			(point.y-otherpoint.y)*(point.y-otherpoint.y)+0.2*R[rank]*R[i]);
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.x-otherpoint.x)/(dist*dist*dist);
		dx+=temp;
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.y-otherpoint.y)/(dist*dist*dist);
		dy+=temp;
		// Right Up side
		otherpoint.y=2*outer_.GetGridBoundary(Up)-otherpoint.y;
		dist=sqrt((point.x-otherpoint.x)*(point.x-otherpoint.x)+
			(point.y-otherpoint.y)*(point.y-otherpoint.y)+0.2*R[rank]*R[i]);
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.x-otherpoint.x)/(dist*dist*dist);
		dx+=temp;
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.y-otherpoint.y)/(dist*dist*dist);
		dy+=temp;
		// Right Down side
		otherpoint.y=2*outer_.GetGridBoundary(Down)-tproc.GetMeshPoint(i).y;
		dist=sqrt((point.x-otherpoint.x)*(point.x-otherpoint.x)+
			(point.y-otherpoint.y)*(point.y-otherpoint.y)+0.2*R[rank]*R[i]);
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.x-otherpoint.x)/(dist*dist*dist);
		dx+=temp;
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.y-otherpoint.y)/(dist*dist*dist);
		dy+=temp;
		// Center bottom side
		otherpoint.x=tproc.GetMeshPoint(i).x;
		dist=sqrt((point.x-otherpoint.x)*(point.x-otherpoint.x)+
			(point.y-otherpoint.y)*(point.y-otherpoint.y)+0.2*R[rank]*R[i]);
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.x-otherpoint.x)/(dist*dist*dist);
		dx+=temp;
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.y-otherpoint.y)/(dist*dist*dist);
		dy+=temp;
		// Left Bottom side
		otherpoint.x=2*outer_.GetGridBoundary(Left)-otherpoint.x;
		dist=sqrt((point.x-otherpoint.x)*(point.x-otherpoint.x)+
			(point.y-otherpoint.y)*(point.y-otherpoint.y)+0.2*R[rank]*R[i]);
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.x-otherpoint.x)/(dist*dist*dist);
		dx+=temp;
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.y-otherpoint.y)/(dist*dist*dist);
		dy+=temp;
		// Left side
		otherpoint.y=tproc.GetMeshPoint(i).y;
		dist=sqrt((point.x-otherpoint.x)*(point.x-otherpoint.x)+
			(point.y-otherpoint.y)*(point.y-otherpoint.y)+0.2*R[rank]*R[i]);
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.x-otherpoint.x)/(dist*dist*dist);
		dx+=temp;
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.y-otherpoint.y)/(dist*dist*dist);
		dy+=temp;
		// Left top side
		otherpoint.y=2*outer_.GetGridBoundary(Up)-otherpoint.y;
		dist=sqrt((point.x-otherpoint.x)*(point.x-otherpoint.x)+
			(point.y-otherpoint.y)*(point.y-otherpoint.y)+0.2*R[rank]*R[i]);
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.x-otherpoint.x)/(dist*dist*dist);
		dx+=temp;
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.y-otherpoint.y)/(dist*dist*dist);
		dy+=temp;
		// Top side
		otherpoint.x=tproc.GetMeshPoint(i).x;
		dist=sqrt((point.x-otherpoint.x)*(point.x-otherpoint.x)+
			(point.y-otherpoint.y)*(point.y-otherpoint.y)+0.2*R[rank]*R[i]);
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.x-otherpoint.x)/(dist*dist*dist);
		dx+=temp;
		temp=M_PI*R[rank]*R[rank]*R[rank]*(PointsPerProc_/NPerProc[i]-1)*
			(point.y-otherpoint.y)/(dist*dist*dist);
		dy+=temp;
	}*/
	// Moving according to pressure
	vector<int> neigh=tproc.GetNeighbors(rank);
	for(int i=0;i<(int)neigh.size();++i)
	{
		if(neigh[i]==-1)
			continue;
		Vector2D otherpoint=tproc.GetMeshPoint(neigh[i]);
		dx+=PointsPerProc_*(1.0/NPerProc[rank]-1.0/NPerProc[neigh[i]])*(
			otherpoint.x-point.x);
		dy+=PointsPerProc_*(1.0/NPerProc[rank]-1.0/NPerProc[neigh[i]])*(
			otherpoint.y-point.y);
	}


	dx=(dx>0) ? min(dx,speed_*R[rank]) : -min(-dx,speed_*R[rank]);
	dy=(dy>0) ? min(dy,speed_*R[rank]) : -min(-dy,speed_*R[rank]);
	dx+=dxround;
	dy+=dyround;
	// Make sure not out of bounds
	const double close=0.99999;
	const double wx=outer_.GetGridBoundary(Right)-outer_.GetGridBoundary(Left);
	const double wy=outer_.GetGridBoundary(Up)-outer_.GetGridBoundary(Down);
	if(point.x+dx>close*outer_.GetGridBoundary(Right))
		if(outer_.GetGridBoundary(Right)-point.x<(1-close)*wx)
			dx=-wx*(1-close);
		else
			dx=0.5*(outer_.GetGridBoundary(Right)-point.x);
	if(point.x+dx<close*outer_.GetGridBoundary(Left))
		if(-outer_.GetGridBoundary(Left)+point.x<(1-close)*wx)
			dx=wx*(1-close);
		else
			dx=-0.5*(-outer_.GetGridBoundary(Left)+point.x);
	if(point.y+dy>close*outer_.GetGridBoundary(Up))
		if(outer_.GetGridBoundary(Up)-point.y<(1-close)*wy)
			dy=-wy*(1-close);
		else
			dy=0.5*(outer_.GetGridBoundary(Up)-point.y);
	if(point.y+dy<close*outer_.GetGridBoundary(Down))
		if(-outer_.GetGridBoundary(Down)+point.y<(1-close)*wy)
			dy=wy*(1-close);
		else
			dy=0.5*(outer_.GetGridBoundary(Down)-point.y);
	vector<Vector2D> cor=tproc.GetMeshPoints();
	cor[rank]=cor[rank]+Vector2D(dx,dy);
	cor.resize(nproc);
	// Have all processors have the same points
	vector<double> dxtemp(nproc),dytemp(nproc);
	MPI_Gather(&cor[rank].x,1,MPI_DOUBLE,&dxtemp[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&dxtemp[0],nproc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&cor[rank].y,1,MPI_DOUBLE,&dytemp[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&dytemp[0],nproc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	for(int i=0;i<nproc;++i)
		cor[i]=Vector2D(dxtemp[i],dytemp[i]);
	tproc.Update(cor);
#endif
}
