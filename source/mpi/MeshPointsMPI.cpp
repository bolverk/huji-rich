#include "MeshPointsMPI.hpp"
#include "../misc/mesh_generator.hpp"
#ifdef RICH_MPI
typedef boost::mt19937_64 gen_type;

namespace
{
	// result it : minx, maxx, miny, maxy
	boost::array<double,4> FindMaxEdges(Tessellation const& tess)
	{
		int rank=0;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		const vector<int> edge_index=tess.GetCellEdges(rank);
		const int n=(int)edge_index.size();
		boost::array<double,4> res;
		res[0]=min(tess.GetEdge(edge_index[0]).vertices.first.x,
			tess.GetEdge(edge_index[0]).vertices.second.x);
		res[1]=max(tess.GetEdge(edge_index[0]).vertices.first.x,
			tess.GetEdge(edge_index[0]).vertices.second.x);
		res[2]=min(tess.GetEdge(edge_index[0]).vertices.first.y,
			tess.GetEdge(edge_index[0]).vertices.second.y);
		res[3]=max(tess.GetEdge(edge_index[0]).vertices.first.y,
			tess.GetEdge(edge_index[0]).vertices.second.y);
		for(int i=1;i<n;++i)
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

	vector<Vector2D> polygon_clip(vector<Vector2D> const& points,
		vector<Vector2D> const& vertices)
	{
		vector<Vector2D> res;
		for(vector<Vector2D>::const_iterator point=points.begin(), end_position=points.end();
			point!=end_position;++point){
				if(PointInCell(vertices,*point))
					res.push_back(*point);
		}
		return res;
	}

	//returns minr,maxr,minangle,maxangle
	boost::array<double,4> GetBindingArc(vector<Vector2D> const& cpoints,
		Vector2D const& center,Vector2D const& procpoint)
	{
		double mincellR=abs(cpoints[0]-center);
		double maxcellR(mincellR);
		int npoints=(int)cpoints.size();
		for(int i=1;i<npoints;++i)
		{
			double temp=abs(cpoints[i]-center);
			maxcellR=max(maxcellR,temp);
			mincellR=min(mincellR,temp);
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
			boost::array<Vector2D,3> tocheck;
			tocheck[0]=center;
			tocheck[1]=procpoint;
			tocheck[2]=cpoints[0];
			double mintemp=orient2d(tocheck)/abs(cpoints[0]-center);
			double maxtemp=mintemp;
			int minloc=0,maxloc=0;
			for(int i=1;i<npoints;++i)
			{
				tocheck[2]=cpoints[i];
				double temp=orient2d(tocheck)/abs(cpoints[i]-center);
				if(temp<mintemp)
				{
					mintemp=temp;
					minloc=i;
				}
				if(temp>maxtemp)
				{
					maxtemp=temp;
					maxloc=i;
				}
			}
			minangle=atan2(cpoints[minloc].y-center.y,cpoints[minloc].x-center.x);
			minangle=(minangle<0) ? (minangle+2*M_PI) : minangle;
			maxangle=atan2(cpoints[maxloc].y-center.y,cpoints[maxloc].x-center.x);
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
	const double Area=(upperright.x-lowerleft.x)*(upperright.y-lowerleft.y);
	const boost::array<double,4> tessEdges=FindMaxEdges(tess);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	const double myarea=tess.GetVolume(rank);
	int mypoints=(int)floor(npoints*myarea/Area+0.5);
	vector<Vector2D> res;
	res.reserve(mypoints);
	vector<Vector2D> cpoints;
	ConvexHull(cpoints,&tess,rank);
	double ran[2];
	gen_type gen(rank);
	boost::random::uniform_real_distribution<> dist;
	// change aboev to have seed==rank
	while((int)res.size()<mypoints)
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
	const double widthx = (upperright.x-lowerleft.x)/(double)nx;
	const double widthy = (upperright.y-lowerleft.y)/(double)ny;
	const boost::array<double,4> tessEdges=FindMaxEdges(tess);
	nx=(int)floor((tessEdges[1]-tessEdges[0])/widthx+0.5);
	ny=(int)floor((tessEdges[3]-tessEdges[2])/widthy+0.5);
	vector<Vector2D> res;
	res.reserve(nx*ny);
	const int nx0=(int)floor((tessEdges[0]-lowerleft.x)/widthx+0.5);
	const int ny0=(int)floor((tessEdges[2]-lowerleft.y)/widthy+0.5);
	Vector2D point;
	int rank=0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	vector<Vector2D> cpoints;
	ConvexHull(cpoints,&tess,rank);
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			point.x = ((double)i+0.5+nx0)*widthx+lowerleft.x;
			point.y = ((double)j+0.5+ny0)*widthy+lowerleft.y;
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
	Vector2D const& bottomleft,Vector2D const& topright,
	Tessellation const& tess,double xc,double yc)
{
	double A=sqrt(M_PI*(Rmax*Rmax-Rmin*Rmin)/PointNum);
	int Nr=int((Rmax-Rmin)/A);
	double dr=(Rmax-Rmin)/Nr;
	int rank=0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	vector<Vector2D> cpoints;
	ConvexHull(cpoints,&tess,rank);
	boost::array<double,4> arc=GetBindingArc(cpoints,Vector2D(xc,yc),
		tess.GetMeshPoint(rank));
	double mincellR=max(min(arc[0],Rmax),Rmin);
	double maxcellR=max(min(arc[1],Rmax),Rmin);
	double minangle=arc[2];
	double maxangle=arc[3];
	int nrmin=(int)((mincellR-Rmin)/dr);
	int nrmax=(int)((maxcellR-Rmin)/dr+0.5);
	vector<Vector2D> res;
	for(int i=nrmin;i<nrmax;++i)
	{
		double r=Rmin+i*dr;
		double dphi=A/r;
		int phimin=(int)(minangle/dphi);
		int phimax=(int)(maxangle/dphi+0.5);
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
	double dphi;
	int Nphi;
	Vector2D pos;
	vector<Vector2D> res;
	int rank=0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	vector<Vector2D> cpoints;
	ConvexHull(cpoints,&tess,rank);
	boost::array<double,4> arc=GetBindingArc(cpoints,Vector2D(xc,yc),
		tess.GetMeshPoint(rank));
	double mincellR=max(min(arc[0],Rmax),Rmin);
	double maxcellR=max(min(arc[1],Rmax),Rmin);
	double minangle=arc[2];
	double maxangle=arc[3];
	double r;
	int nrmin=(int)((pow(max(mincellR,Rmin),alpha+1)-pow(Rmin,alpha+1))/(2*M_PI*(alpha+1)/N0));
	int nrmax=(int)((pow(min(maxcellR,Rmax),alpha+1)-pow(Rmin,alpha+1))/(2*M_PI*(alpha+1)/N0)+0.5);
	for(int i=nrmin;i<nrmax;++i)
	{
		r=pow(2*M_PI*i*(alpha+1)/N0+pow(Rmin,alpha+1),1.0/(alpha+1));
		Nphi=int(floor(N0*pow(r,1+alpha)+1.5));
		dphi=2*M_PI/Nphi;
		int phimin=(int)(minangle/dphi);
		int phimax=(int)(maxangle/dphi+0.5);
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
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	vector<Vector2D> cpoints;
	ConvexHull(cpoints,&tproc,rank);
	for(int i=0;i<point_number;++i)
	{
		const double angle = 2*M_PI*double(i)/double(point_number);
		Vector2D temp= center+pol2cart(radius,angle);
		if(PointInCell(cpoints,temp))
			res.push_back(temp);
	}
	return res;
}


namespace {
	vector<size_t> most_uniform_range_partition(size_t range,
		size_t partition_number)
	{
		vector<size_t> res(partition_number, range/partition_number);
		for(size_t i=partition_number-range%partition_number;
			i<partition_number;++i)
			res[i] += 1;
		return res;
	}

	size_t sort_point(Tessellation const& tess,
		Vector2D const& point)
	{
		for(size_t i=0;i<tess.GetPointNo();++i){
			vector<Vector2D> vertices;
			ConvexHull(vertices,&tess,i);
			if(PointInCell(vertices,point))
				return i;
		}
		assert(false && "Point does not belong in tessllation");
	}

	vector<vector<Vector2D> > sort_points
		(Tessellation const& process_tess,
		Index2Member<Vector2D> const& grid_generator,
		size_t start, size_t ending)
	{
		vector<vector<Vector2D> > res(process_tess.GetPointNo());
		for(size_t i=start;i<ending;++i)
		{
			const Vector2D point = grid_generator(i);
			res[sort_point(process_tess,point)].push_back(point);
		}
		return res;
	}
}

vector<Vector2D> distribute_grid(Tessellation const& process_tess,
	Index2Member<Vector2D> const& grid_generator)
{
	assert(get_mpi_size()==process_tess.GetPointNo() &&
		"Number of processors must be equal to the number of cells");
	const vector<size_t> range_list = most_uniform_range_partition
		(grid_generator.getLength(), process_tess.GetPointNo());
	size_t start = 0;
	for(size_t i=0;i<get_mpi_rank();++i)
		start += range_list[i];
	size_t ending = start + range_list[get_mpi_rank()];
	const vector<vector<Vector2D> > sorted_points = 
		sort_points(process_tess, grid_generator, start, ending);
	vector<Vector2D> res = sorted_points[get_mpi_rank()];
	for(size_t i=0;i<get_mpi_rank();++i){
		MPI_VectorSend_Vector2D(sorted_points[i],i,0,MPI_COMM_WORLD);
		vector<Vector2D> buf;
		MPI_VectorRecv_Vector2D(buf,i,0,MPI_COMM_WORLD);
		res.reserve(res.size()+distance(buf.begin(),buf.end()));
		res.insert(res.end(),buf.begin(),buf.end());
	}
	for(size_t i=get_mpi_rank()+1;i<get_mpi_size();++i){
		vector<Vector2D> buf;
		MPI_VectorRecv_Vector2D(buf,i,0,MPI_COMM_WORLD);
		MPI_VectorSend_Vector2D(sorted_points[i],i,0,MPI_COMM_WORLD);
		res.reserve(res.size()+distance(buf.begin(),buf.end()));
		res.insert(res.end(),buf.begin(),buf.end());
	}
	return res;
}

CartesianGridGenerator::CartesianGridGenerator
	(size_t nx, size_t ny, const Vector2D& lower_left, const Vector2D& upper_right):
nx_(nx), ny_(ny), lower_left_(lower_left), upper_right_(upper_right) {}

size_t CartesianGridGenerator::getLength(void) const
{
	return nx_*ny_;
}

Vector2D CartesianGridGenerator::operator()(size_t idx) const
{
	const size_t i = idx/nx_;
	const size_t j = idx%ny_;
	return lower_left_ +
		Vector2D((upper_right_-lower_left_).x*(0.5+(double)i)/(double)nx_,
		(upper_right_-lower_left_).y*(0.5+(double)j)/(double)ny_);
}

#endif
