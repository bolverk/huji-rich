#include "RoundGrid.hpp"
#include "HilbertOrder.hpp"

vector<Vector2D> RoundGrid(vector<Vector2D> const& points,
			   OuterBoundary const* bc,int NumberIt,
			   #ifdef RICH_MPI
			   Tessellation const* tproc,
			   #endif
	Tessellation *tess)
{
	vector<int> indeces = HilbertOrder(points, static_cast<int>(points.size()));
	vector<Vector2D> res(points);
	res = VectorValues(res, indeces);
	VoronoiMesh default_tess;
	if(tess==0)
		tess=&default_tess;
#ifdef RICH_MPI
	if(tproc==0)
		tess->Initialise(points,bc);
	else
		tess->Initialise(points,*tproc,bc);
#else
	tess->Initialise(points,bc);
#endif
	double pi= 3.141592653;
	double eta_=0.02,chi_=1;
	int N=tess->GetPointNo();

	// Copy the points
	for(int i=0;i<N;++i)
	  res[static_cast<size_t>(i)]=tess->GetMeshPoint(i);

	for(int j=0;j<NumberIt;++j)
	{
#ifdef RICH_MPI
		N=tess->GetPointNo();
		res=tess->GetMeshPoints();
		res.resize(static_cast<size_t>(N));
#endif
		for(int i=0;i<N;++i)
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
			res[static_cast<size_t>(i)]=tess->GetMeshPoint(i)+dw;
		}
#ifdef RICH_MPI
		if(tproc==0)
			tess->Update(res);
		else
			tess->Update(res,*tproc);
#else
		tess->Update(res);
#endif
	}
#ifdef RICH_MPI
	N=tess->GetPointNo();
	res=tess->GetMeshPoints();
	res.resize(static_cast<size_t>(N));
#endif
	return res;
}
