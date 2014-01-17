#include "RoundGrid.hpp"
vector<Vector2D> RoundGrid(vector<Vector2D> const& points,
	OuterBoundary *bc,int NumberIt,int InnerNum,Tessellation *tess)
{
	VoronoiMesh default_tess;
	if(tess==0)
		tess=&default_tess;
	tess->Initialise(points,bc);
	double pi= 3.141592653;
	double eta_=0.05,chi_=1;
	int N=tess->GetPointNo();

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
			res[i]=tess->GetMeshPoint(i)+dw;
		}
		tess->Update(res);
	}
	
	return res;
}
