#include "RoundGrid.hpp"
vector<Vector2D> RoundGrid(vector<Vector2D> const& points,
	OuterBoundary const* bc,int NumberIt,int InnerNum,Tessellation const* tproc,
	Tessellation *tess)
{
	VoronoiMesh default_tess;
	if(tess==0)
		tess=&default_tess;
	tess->Initialise(points,bc);
	double pi= 3.141592653;
	double eta_=0.05,chi_=1;
	int N=tess->GetPointNo();
	vector<Vector2D> cpoints;
#ifdef RICH_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	ConvexHull(cpoints,tproc,rank);
#endif
	// Copy the points
	vector<Vector2D> res(N);
	for(int i=0;i<N;++i)
		res[i]=tess->GetMeshPoint(i);

	for(int j=0;j<NumberIt;++j)
	{
		for(int i=InnerNum;i<N;++i)
		{
			double R = sqrt(tess->GetVolume(i)/pi);
			Vector2D s = tess->GetCellCM(i);
			Vector2D r = tess->GetMeshPoint(i);
			double d = abs(s-r);
			Vector2D dw;
			if(d/eta_/R<0.95)
				dw = 0*s;
			else 
				dw = chi_*0.5*(s-r);
#ifdef RICH_MPI
			if(!PointInCell(cpoints,tess->GetMeshPoint(i)+dw))
				dw=Vector2D(0,0);
#endif
			res[i]=tess->GetMeshPoint(i)+dw;
		}
		tess->Update(res);
	}
	
	return res;
}
