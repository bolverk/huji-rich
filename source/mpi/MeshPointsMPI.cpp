#include "MeshPointsMPI.hpp"
#include "../misc/mesh_generator.hpp"
#ifdef RICH_MPI
#include <mpi.h>
typedef boost::mt19937_64 gen_type;

namespace
{
	// result it : minx, maxx, miny, maxy
	boost::array<double,4> FindMaxEdges(Tessellation const& tess)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		const vector<int>& edge_index=tess.GetCellEdges(rank);
		const int n=static_cast<int>(edge_index.size());
		boost::array<double,4> res;
		res[0]=min(tess.GetEdge(edge_index[0]).vertices.first.x,
			tess.GetEdge(edge_index[0]).vertices.second.x);
		res[1]=max(tess.GetEdge(edge_index[0]).vertices.first.x,
			tess.GetEdge(edge_index[0]).vertices.second.x);
		res[2]=min(tess.GetEdge(edge_index[0]).vertices.first.y,
			tess.GetEdge(edge_index[0]).vertices.second.y);
		res[3]=max(tess.GetEdge(edge_index[0]).vertices.first.y,
			tess.GetEdge(edge_index[0]).vertices.second.y);
		for(size_t i=1;i<static_cast<size_t>(n);++i)
		{
			res[0]=min(min(tess.GetEdge(edge_index[i]).vertices.first.x,
				tess.GetEdge(edge_index[i]).vertices.second.x),res[0]);
			res[1]=max(max(tess.GetEdge(edge_index[i]).vertices.first.x,
				tess.GetEdge(edge_index[i]).vertices.second.x),res[1]);
			res[2]=min(min(tess.GetEdge(edge_index[i]).vertices.first.y,
				tess.GetEdge(edge_index[i]).vertices.second.y),res[2]);
			res[3]=max(max(tess.GetEdge(edge_index[i]).vertices.first.y,
				tess.GetEdge(edge_index[i]).vertices.second.y),res[3]);
		}
		return res;
	}

	//returns minr,maxr,minangle,maxangle
	boost::array<double,4> GetBindingArc(vector<Vector2D> const& cpoints,
		Vector2D const& center,Vector2D const& procpoint)
	{
		Edge etemp;
		etemp.vertices.first=cpoints[0];
		etemp.vertices.second=cpoints[1];
		double mincellR=abs(cpoints[0]-center);
		double maxcellR(mincellR);
		mincellR=min(mincellR,DistanceToEdge(center,etemp));
		int npoints=static_cast<int>(cpoints.size());
		for(size_t i=1;i<static_cast<size_t>(npoints);++i)
		{
			double temp=abs(cpoints[i]-center);
			maxcellR=max(maxcellR,temp);
			etemp.vertices.first=cpoints[i];
			etemp.vertices.second=cpoints[(i+1)%static_cast<size_t>(npoints)];
			mincellR=min(min(mincellR,temp),DistanceToEdge(center,etemp));
		}
		// Get the maximum and minimum angles
		double maxangle,minangle;
		if(PointInCell(cpoints,center))
		{
			mincellR=0;
			maxangle=2*M_PI;
			minangle=0;
		}
		else
		{
			const Vector2D tocenter=procpoint-center;
			const double Rcenter=abs(tocenter);
			boost::array<Vector2D,3> tocheck;
			tocheck[0]=center;
			tocheck[1]=procpoint;
			double minusval=0,plusval=0;
			int minusloc=0,plusloc=0;
			for(size_t i=0;i<static_cast<size_t>(npoints);++i)
			{
				tocheck[2]=cpoints[i];
				const Vector2D tocpoint=cpoints[i]-center;
				const double angle=acos(ScalarProd(tocenter,tocpoint)/(abs(tocpoint)*Rcenter));
				if(orient2d(TripleConstRef<Vector2D>(center,
								     procpoint,
								     cpoints[i]))>0)
				{	
					if(angle>plusval)
					{
					  plusloc=static_cast<int>(i);
						plusval=angle;
					}
				}
				else
				{
					if(angle>minusval)
					{
					  minusloc=static_cast<int>(i);
						minusval=angle;
					}
				}
			}
			minangle=atan2(cpoints[static_cast<size_t>(minusloc)].y-center.y,cpoints[static_cast<size_t>(minusloc)].x-center.x);
			minangle=(minangle<0) ? (minangle+2*M_PI) : minangle;
			maxangle=atan2(cpoints[static_cast<size_t>(plusloc)].y-center.y,cpoints[static_cast<size_t>(plusloc)].x-center.x);
			maxangle=(maxangle<0) ? (maxangle+2*M_PI) : maxangle;
			if(maxangle<minangle)
				minangle-=2*M_PI;
		}
		boost::array<double,4> res;
		res[0]=mincellR;
		res[1]=maxcellR;
		res[2]=minangle;
		res[3]=maxangle;
		return res;
	}
}

vector<Vector2D> RandSquare(int npoints,Tessellation const& tess,
	Vector2D const& lowerleft,Vector2D const& upperright)
{
	boost::array<double,4> tessEdges=FindMaxEdges(tess);
	if (tessEdges[0]>upperright.x || tessEdges[1]<lowerleft.x || tessEdges[2]>upperright.y ||
		tessEdges[3] < lowerleft.y)
		return vector<Vector2D>();
	tessEdges[0] = std::max(tessEdges[0], lowerleft.x);
	tessEdges[2] = std::max(tessEdges[2], lowerleft.y);
	tessEdges[1] = std::min(tessEdges[1], upperright.x);
	tessEdges[3] = std::min(tessEdges[3], upperright.y);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const double myarea=std::min(tess.GetVolume(rank), (tessEdges[1]- tessEdges[0])*(tessEdges[3] - tessEdges[2]));
	const double Area = (upperright.x - lowerleft.x)*(upperright.y - lowerleft.y);
	int mypoints=static_cast<int>(floor(npoints*myarea/Area+0.5));
	vector<Vector2D> res;
	res.reserve(static_cast<size_t>(mypoints));
	vector<Vector2D> cpoints;
	ConvexHull(cpoints,tess,rank);
	double ran[2];
	gen_type gen(static_cast<size_t>(rank));
	boost::random::uniform_real_distribution<> dist;
	// change aboev to have seed==rank
	while(static_cast<int>(res.size())<mypoints)
	{
		ran[0]=dist(gen)*(tessEdges[1]-tessEdges[0])+tessEdges[0];
		ran[1]=dist(gen)*(tessEdges[3]-tessEdges[2])+tessEdges[2];
		Vector2D p(ran[0],ran[1]);
		if(PointInCell(cpoints,p))
			res.push_back(p);
	}
	return res;
}

vector<Vector2D> SquareMeshM(int nx,int ny,Tessellation const& tess,
	Vector2D const&lowerleft,Vector2D const&upperright)
{
	const double widthx = (upperright.x-lowerleft.x)/static_cast<double>(nx);
	const double widthy = (upperright.y-lowerleft.y)/static_cast<double>(ny);
	const boost::array<double,4> tessEdges=FindMaxEdges(tess);
	nx=static_cast<int>(floor((tessEdges[1]-tessEdges[0])/widthx+2));
	ny=static_cast<int>(floor((tessEdges[3]-tessEdges[2])/widthy+2));
	vector<Vector2D> res;
	res.reserve(static_cast<size_t>(nx*ny));
	const int nx0=static_cast<int>(floor((tessEdges[0]-lowerleft.x)/widthx-0.5));
	const int ny0=static_cast<int>(floor((tessEdges[2]-lowerleft.y)/widthy-0.5));
	Vector2D point;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<Vector2D> cpoints;
	ConvexHull(cpoints,tess,rank);
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
		  point.x = (static_cast<double>(i)+0.5+nx0)*widthx+lowerleft.x;
		  point.y = (static_cast<double>(j)+0.5+ny0)*widthy+lowerleft.y;
			if((point.x<lowerleft.x)||(point.x>upperright.x)||
				(point.y<lowerleft.y)||(point.y>upperright.y))
				continue;
			if(PointInCell(cpoints,point))
				res.push_back(point);
		}
	}
	return res;
}

vector<Vector2D> CirclePointsRmaxM(int PointNum,double Rmin,double Rmax,
				   Vector2D const& /*bottomleft*/,
				   Vector2D const& /*topright*/,
	Tessellation const& tess,double xc,double yc)
{
	double A=sqrt(M_PI*(Rmax*Rmax-Rmin*Rmin)/PointNum);
	int Nr=int((Rmax-Rmin)/A);
	double dr=(Rmax-Rmin)/Nr;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<Vector2D> cpoints;
	ConvexHull(cpoints,tess,rank);
	boost::array<double,4> arc=GetBindingArc(cpoints,Vector2D(xc,yc),
		tess.GetMeshPoint(rank));
	double mincellR=max(min(arc[0],Rmax),Rmin);
	double maxcellR=max(min(arc[1],Rmax),Rmin);
	double minangle=arc[2];
	double maxangle=arc[3];
	int nrmin=static_cast<int>((mincellR-Rmin)/dr);
	int nrmax=static_cast<int>((maxcellR-Rmin)/dr+0.5);
	vector<Vector2D> res;
	for(int i=nrmin;i<nrmax;++i)
	{
		double r=Rmin+i*dr;
		double dphi=A/r;
		int phimin=static_cast<int>(minangle/dphi);
		int phimax=static_cast<int>(maxangle/dphi+0.5);
		for(int j=phimin;j<phimax;++j)
		{
			Vector2D temp(Vector2D(r*cos(dphi*j)+xc,r*sin(dphi*j)+yc));
			if(PointInCell(cpoints,temp))
				res.push_back(temp);
		}
	}
	return res;
}

vector<Vector2D> CirclePointsRmax_aM(int PointNum,double Rmin,double Rmax,
	double xc,double yc,double alpha,Tessellation const& tess)
{
	double N0=sqrt(PointNum*4*M_PI*(alpha+1)/(pow(Rmax,2*(alpha+1))-
		pow(Rmin,2*(alpha+1))));
	Vector2D pos;
	vector<Vector2D> res;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<Vector2D> cpoints;
	ConvexHull(cpoints,tess,rank);
	boost::array<double,4> arc=GetBindingArc(cpoints,Vector2D(xc,yc),
		tess.GetMeshPoint(rank));
	double mincellR=max(min(arc[0],Rmax),Rmin);
	double maxcellR=max(min(arc[1],Rmax),Rmin);
	double minangle=arc[2];
	double maxangle=arc[3];
	int nrmin=static_cast<int>((pow(max(mincellR,Rmin),alpha+1)-pow(Rmin,alpha+1))/(2*M_PI*(alpha+1)/N0));
	int nrmax=static_cast<int>((pow(min(maxcellR,Rmax),alpha+1)-pow(Rmin,alpha+1))/(2*M_PI*(alpha+1)/N0)+0.5);
	for(int i=nrmin;i<nrmax;++i)
	{
		const double r=pow(2*M_PI*i*(alpha+1)/N0+pow(Rmin,alpha+1),1.0/(alpha+1));
		const int Nphi=int(floor(N0*pow(r,1+alpha)+1.5));
		const double dphi=2*M_PI/Nphi;
		int phimin=static_cast<int>(minangle/dphi);
		int phimax=static_cast<int>(maxangle/dphi+0.5);
		for(int j=phimin;j<phimax;++j)
		{
			pos.Set(r*cos(dphi*j)+xc,r*sin(dphi*j)+yc);
			if(PointInCell(cpoints,pos))
				res.push_back(pos);
		}
	}
	return res;
}

vector<Vector2D> circle_circumferenceM(int point_number,double radius,
	Vector2D const& center,Tessellation const& tproc)
{
	vector<Vector2D> res;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<Vector2D> cpoints;
	ConvexHull(cpoints,tproc,rank);
	for(int i=0;i<point_number;++i)
	{
		const double angle = 2*M_PI*double(i)/double(point_number);
		Vector2D temp= center+pol2cart(radius,angle);
		if(PointInCell(cpoints,temp))
			res.push_back(temp);
	}
	return res;
}

#endif
