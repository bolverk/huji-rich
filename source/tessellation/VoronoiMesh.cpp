#include "VoronoiMesh.hpp"
#include <cmath>
#include "../misc/simple_io.hpp"
#include "hdf5_logger.hpp"

using std::abs;

namespace 
{
#ifdef RICH_MPI

  template<class T> void tidy(vector<T>& v)
  {
    if(!v.empty()){
      sort(v.begin(),v.end());
      v = unique(v);
    }
  }

  vector<Vector2D> my_convex_hull(const Tessellation& tess,
				  int index)
  {
    vector<Vector2D> res;
    ConvexHull(res,&tess,index);
    return res;
  }

	// Send/Recv points from other processors and moves periodic points.
	// Returns the new points as well as the self index of the points that where kept
	// and the indeces of the sent points
	vector<Vector2D> UpdateMPIPoints(Tessellation const& vproc,int rank,
		vector<int> cornerproc,vector<int> const& proclist,vector<Vector2D>
		const& points,OuterBoundary const* obc,vector<int> &selfindex,
		vector<int> &sentproc,vector<vector<int> > &sentpoints)
	{
	  tidy(cornerproc);
		vector<Vector2D> res;
		res.reserve(points.size());
		selfindex.clear();
		int npoints=static_cast<int>(points.size());
		int nproc=vproc.GetPointNo();
		const double dx=obc->GetGridBoundary(Right)-obc->GetGridBoundary(Left);
		const double dy=obc->GetGridBoundary(Up)-obc->GetGridBoundary(Down);
		const vector<Vector2D> cproc = my_convex_hull(vproc,rank);
		int ncorner=static_cast<int>(cornerproc.size());
		vector<vector<Vector2D> > cornerpoints;
		vector<int> realcornerproc;
		for(int i=0;i<ncorner;++i)
		{
			if(cornerproc[static_cast<size_t>(i)]>-1)
			{
				realcornerproc.push_back(vproc.GetOriginalIndex(cornerproc[static_cast<size_t>(i)]));
				cornerpoints.push_back(vector<Vector2D> ());
				ConvexHull(cornerpoints.back(),&vproc,realcornerproc.back());
			}
		}
		ncorner=static_cast<int>(realcornerproc.size());
		vector<vector<Vector2D> > cornersend(cornerpoints.size());
		int nneigh=static_cast<int>(proclist.size());
		vector<vector<Vector2D> > neighpoints;
		vector<int> realneighproc;
		for(int i=0;i<nneigh;++i)
		{
			if(proclist[static_cast<size_t>(i)]>-1)
			{
				realneighproc.push_back(vproc.GetOriginalIndex(proclist[static_cast<size_t>(i)]));
				neighpoints.push_back(vector<Vector2D> ());
				ConvexHull(neighpoints.back(),&vproc,realneighproc.back());
			}
		}
		nneigh=static_cast<int>(realneighproc.size());
		sentpoints.clear();
		sentproc=realneighproc;
		sentproc.insert(sentproc.end(),realcornerproc.begin(),realcornerproc.end());
		sentpoints.resize(sentproc.size());
		vector<vector<Vector2D> > neighsend(realneighproc.size());
		int pointoutside=0;
		for(int i=0;i<npoints;++i)
		{
			Vector2D temp(points[static_cast<size_t>(i)]);
			if(obc->GetBoundaryType()==Periodic)
			{
				if(temp.x>obc->GetGridBoundary(Right))
					temp.x-=dx;
				if(temp.x<obc->GetGridBoundary(Left))
					temp.x+=dx;
				if(temp.y>obc->GetGridBoundary(Up))
					temp.y-=dy;
				if(temp.y<obc->GetGridBoundary(Down))
					temp.y+=dy;
			}
			if(obc->GetBoundaryType()==HalfPeriodic)
			{
				if(temp.x>obc->GetGridBoundary(Right))
					temp.x-=dx;
				if(temp.x<obc->GetGridBoundary(Left))
					temp.x+=dx;
			}
			if(PointInCell(cproc,temp))
			{
				res.push_back(temp);
				selfindex.push_back(i);
				continue;
			}
			bool good=false;
			for(int j=0;j<nneigh;++j)
			{
				if(PointInCell(neighpoints[static_cast<size_t>(j)],temp))
				{
					neighsend[static_cast<size_t>(j)].push_back(temp);
					sentpoints[static_cast<size_t>(j)].push_back(i);
					good=true;
					break;
				}
			}
			if(good)
				continue;
			for(int j=0;j<ncorner;++j)
			{
				if(PointInCell(cornerpoints[static_cast<size_t>(j)],temp))
				{
					cornersend[static_cast<size_t>(j)].push_back(temp);
					sentpoints[static_cast<size_t>(j+nneigh)].push_back(i);
					good=true;
					break;
				}
			}
			if(good)
				continue;
			for(int j=0;j<nproc;++j)
			{
				vector<Vector2D> cellpoints;
				ConvexHull(cellpoints,&vproc,j);
				if(PointInCell(cellpoints,temp))
				{
					bool notadded=false;
					int indexadd=0;
					if(static_cast<int>(realneighproc.size())==nneigh)
						notadded=true;
					if(!notadded)
					{
					  indexadd=static_cast<int>(find(realneighproc.begin()+nneigh,
									 realneighproc.end(),j)-realneighproc.begin());
						if(indexadd==static_cast<int>(realneighproc.size()))
							notadded=true;
					}
					if(notadded)
					{
						realneighproc.insert(realneighproc.end(),j);
						vector<int> vtemp;
						vtemp.push_back(i);
						sentpoints.insert(sentpoints.end(),vtemp);
						vector<Vector2D> v2temp;
						v2temp.push_back(temp);
						neighsend.insert(neighsend.end(),v2temp);
						good=true;
						pointoutside=1;
						break;
					}
					else
					{
						vector<int> vtemp;
						vtemp.push_back(i);
						if(indexadd<nneigh)
							sentpoints[static_cast<size_t>(indexadd)].insert(sentpoints[static_cast<size_t>(indexadd)].end(),vtemp.begin(),
							vtemp.end());
						else
						  sentpoints[static_cast<size_t>(indexadd+ncorner)].insert(sentpoints[static_cast<size_t>(indexadd+ncorner)].end(),vtemp.begin(),
							vtemp.end());
						vector<Vector2D> v2temp;
						v2temp.push_back(temp);
						neighsend[static_cast<size_t>(indexadd)].insert(neighsend[static_cast<size_t>(indexadd)].end(),v2temp.begin(),
							v2temp.end());
						good=true;
						pointoutside=1;
						break;
					}
				}
			}
			if(good)
				continue;
			HDF5Logger log("verror"+int2str(get_mpi_rank())+".h5");
			log.output(vproc);
			UniversalError eo("Point is not inside any processor");
			eo.AddEntry("CPU rank",rank);
			eo.AddEntry("Point number",i);
			eo.AddEntry("Point x cor",points[static_cast<size_t>(i)].x);
			eo.AddEntry("Point y cor",points[static_cast<size_t>(i)].y);
			for(int k=0;k<nneigh;++k){
				eo.AddEntry("Neighbor "+int2str(k)+" is cpu "
					    ,static_cast<double>(realneighproc[static_cast<size_t>(k)]));
				for(size_t l=0;l<neighpoints[static_cast<size_t>(k)].size();++l){
				  eo.AddEntry("neighbor point x",neighpoints[static_cast<size_t>(k)][l].x);
				  eo.AddEntry("neighbor point y",neighpoints[static_cast<size_t>(k)][l].y);
				}
			}			
			for(int k=0;k<ncorner;++k)
				eo.AddEntry("Corner "+int2str(k)+" is cpu "
					    ,static_cast<double>(realcornerproc[static_cast<size_t>(k)]));
			throw eo;
		}
		// Send/Recv the points
		vector<int> procorder=GetProcOrder(rank,nproc);
		// Combine the vectors
		int ntemp=static_cast<int>(realneighproc.size());
		vector<int> neightemp;
		vector<vector<int> > senttemp;
		vector<vector<Vector2D> > tosend;
		int pointoutside2;
		MPI_Allreduce(&pointoutside,&pointoutside2,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
		if(pointoutside2==1)
		{
			vector<int> allsend=join(realneighproc,realcornerproc);
			vector<int> sendnumber(nproc,0),scounts(nproc,1);
			int nsend=static_cast<int>(allsend.size());
			for(int jj=0;jj<nsend;++jj)
			  sendnumber[static_cast<size_t>(allsend[static_cast<size_t>(jj)])]+=1;
			int nrecv;
			MPI_Reduce_scatter(&sendnumber[0],&nrecv,&scounts[0],MPI_INT,MPI_SUM,
				MPI_COMM_WORLD);
			MPI_Request *request=new MPI_Request[static_cast<size_t>(nsend)];
			char temp;
			for(int jj=0;jj<nsend;++jj)
				MPI_Isend(&temp,1,MPI_CHAR,allsend[static_cast<size_t>(jj)],1,MPI_COMM_WORLD,
				&request[static_cast<size_t>(jj)]);
			MPI_Status status;
			vector<int> newtalk;
			for(int jj=0;jj<nrecv;++jj)
			{
				MPI_Probe(MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
				MPI_Recv(&temp,1,MPI_CHAR,status.MPI_SOURCE,1,MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);
				if(find(allsend.begin(),allsend.end(),status.MPI_SOURCE)==
					allsend.end())
					newtalk.push_back(status.MPI_SOURCE);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			delete [] request;
			for(int i=0;i<nneigh;++i)
			{
				neightemp.push_back(realneighproc[static_cast<size_t>(i)]);
				senttemp.push_back(sentpoints[static_cast<size_t>(i)]);
				tosend.push_back(neighsend[static_cast<size_t>(i)]);
			}
			ntemp=static_cast<int>(realcornerproc.size());
			for(int i=0;i<ntemp;++i)
			{
			  int index=static_cast<int>(find(realneighproc.begin(),realneighproc.end(),realcornerproc[static_cast<size_t>(i)])
						     -realneighproc.begin());
				if(index<static_cast<int>(realneighproc.size()))
				{
				  senttemp[static_cast<size_t>(index)].insert(senttemp[static_cast<size_t>(index)].begin(),sentpoints[static_cast<size_t>(i+nneigh)].begin(),
									      sentpoints[static_cast<size_t>(i+nneigh)].end());
					tosend[static_cast<size_t>(index)].insert(tosend[static_cast<size_t>(index)].end(),cornersend[static_cast<size_t>(i)].begin(),
						cornersend[static_cast<size_t>(i)].end());
				}
				else
				{
					neightemp.push_back(realcornerproc[static_cast<size_t>(i)]);
					senttemp.push_back(sentpoints[static_cast<size_t>(i+nneigh)]);
					tosend.push_back(cornersend[static_cast<size_t>(i)]);
				}
			}
			for(int i=nneigh;i<static_cast<int>(realneighproc.size());++i)
			{
				neightemp.push_back(realneighproc[static_cast<size_t>(i)]);
				senttemp.push_back(sentpoints[static_cast<size_t>(i+ntemp)]);
				tosend.push_back(neighsend[static_cast<size_t>(i)]);
			}
			ntemp=static_cast<int>(newtalk.size());
			for(int i=0;i<ntemp;++i)
			{
				senttemp.push_back(vector<int> ());
				neightemp.push_back(newtalk[static_cast<size_t>(i)]);
				tosend.push_back(vector<Vector2D> ());
			}
		}
		else
		{
			for(int i=0;i<ntemp;++i)
			{
				neightemp.push_back(realneighproc[static_cast<size_t>(i)]);
				senttemp.push_back(sentpoints[static_cast<size_t>(i)]);
				tosend.push_back(neighsend[static_cast<size_t>(i)]);
			}
			ntemp=static_cast<int>(realcornerproc.size());
			for(int i=0;i<ntemp;++i)
			{
			  int index=static_cast<int>(find(realneighproc.begin(),realneighproc.end(),realcornerproc[static_cast<size_t>(i)])
						     -realneighproc.begin());
				if(index<static_cast<int>(realneighproc.size()))
				{
				  senttemp[static_cast<size_t>(index)].insert(senttemp[static_cast<size_t>(index)].begin(),sentpoints[static_cast<size_t>(i+nneigh)].begin(),
									      sentpoints[static_cast<size_t>(i+nneigh)].end());
					tosend[static_cast<size_t>(index)].insert(tosend[static_cast<size_t>(index)].end(),cornersend[static_cast<size_t>(i)].begin(),
						cornersend[static_cast<size_t>(i)].end());
				}
				else
				{
					neightemp.push_back(realcornerproc[static_cast<size_t>(i)]);
					senttemp.push_back(sentpoints[static_cast<size_t>(i+nneigh)]);
					tosend.push_back(cornersend[static_cast<size_t>(i)]);
				}
			}
		}
		sentpoints=senttemp;
		sentproc=neightemp;
		vector<Vector2D> toadd=MPI_MassSendRecvVectorVector2D(tosend,sentproc,procorder);
		res.insert(res.end(),toadd.begin(),toadd.end());
		return res;
	}
#endif

	vector<Vector2D> UpdatePoints(vector<Vector2D> const& points,OuterBoundary const* obc)
	{
		if(obc->GetBoundaryType()==Rectengular)
			return points;
		vector<Vector2D> res;
		res.reserve(points.size());
		int npoints=static_cast<int>(points.size());
		const double dx=obc->GetGridBoundary(Right)-obc->GetGridBoundary(Left);
		const double dy=obc->GetGridBoundary(Up)-obc->GetGridBoundary(Down);
		for(int i=0;i<npoints;++i)
		{
			Vector2D temp(points[static_cast<size_t>(i)]);
			if(obc->GetBoundaryType()==Periodic)
			{
				if(temp.x>obc->GetGridBoundary(Right))
					temp.x-=dx;
				if(temp.x<obc->GetGridBoundary(Left))
					temp.x+=dx;
				if(temp.y>obc->GetGridBoundary(Up))
					temp.y-=dy;
				if(temp.y<obc->GetGridBoundary(Down))
					temp.y+=dy;
			}
			if(obc->GetBoundaryType()==HalfPeriodic)
			{
				if(temp.x>obc->GetGridBoundary(Right))
					temp.x-=dx;
				if(temp.x<obc->GetGridBoundary(Left))
					temp.x+=dx;
			}
			res.push_back(temp);
		}
		return res;
	}

#ifdef RICH_MPI
		  /*
	void RemoveDuplicateCorners(vector<int> &cornerproc,vector<vector<int> >
		&corners,vector<int> const& proclist,vector<vector<int> > &toduplicate)
	{
		vector<int> newcornerproc;
		vector<vector<int> > newcorners;
		int ncorners=static_cast<int>(cornerproc.size());
		int nproc=static_cast<int>(proclist.size());
		for(int i=0;i<ncorners;++i)
		{
			if(cornerproc[static_cast<size_t>(i)]==-1)
				continue;
			int index=static_cast<int>(find(proclist.begin(),proclist.end(),cornerproc[static_cast<size_t>(i)])-proclist.begin());
			if(index<nproc)
			{
				// Add unduplicated data
				vector<int> toadd;
				int npoints=static_cast<int>(corners[static_cast<size_t>(i)].size());
				for(int j=0;j<npoints;++j)
				{
					if(!binary_search(toduplicate[static_cast<size_t>(index)].begin(),toduplicate[static_cast<size_t>(index)].end(),
						corners[static_cast<size_t>(i)][static_cast<size_t>(j)]))
						toadd.push_back(corners[static_cast<size_t>(i)][static_cast<size_t>(j)]);
				}
				if(!toadd.empty())
				{
					toduplicate[static_cast<size_t>(index)].insert(toduplicate[static_cast<size_t>(index)].end(),toadd.begin(),
						toadd.end());
					sort(toduplicate[static_cast<size_t>(index)].begin(),toduplicate[static_cast<size_t>(index)].end());
				}
			}
			else
			{
				// copy the data
				newcornerproc.push_back(cornerproc[static_cast<size_t>(i)]);
				newcorners.push_back(corners[static_cast<size_t>(i)]);
			}
		}
		// Remove same CPUs
		if(newcornerproc.size()>1)
		{
			int n=static_cast<int>(newcornerproc.size());
			vector<int> index;
			sort_index(newcornerproc,index);
			sort(newcornerproc.begin(),newcornerproc.end());
			vector<int> temp=unique(newcornerproc);
			int nuinq=static_cast<int>(temp.size());
			vector<vector<int> > cornerstemp(temp.size());
			for(int i=0;i<n;++i)
			{
			  int place=static_cast<int>(find(temp.begin(),temp.end(),newcornerproc[static_cast<size_t>(i)])-temp.begin());
				cornerstemp[static_cast<size_t>(place)].insert(cornerstemp[static_cast<size_t>(place)].begin(),
									       newcorners[static_cast<size_t>(index[static_cast<size_t>(i)])].begin(),newcorners[static_cast<size_t>(index[static_cast<size_t>(i)])].end());
			}
			for(int i=0;i<nuinq;++i)
			{
				sort(cornerstemp[static_cast<size_t>(i)].begin(),cornerstemp[static_cast<size_t>(i)].end());
				cornerstemp[static_cast<size_t>(i)]=unique(cornerstemp[static_cast<size_t>(i)]);
			}
			newcornerproc=temp;
			newcorners=cornerstemp;
		}
		cornerproc=newcornerproc;
		corners=newcorners;
	}
		  */
#endif // RICH_MPI

	Vector2D GetPeriodicDiff(Edge const& edge,OuterBoundary const* obc)
	{
		Vector2D diff;
		if(abs(edge.vertices.first.x-edge.vertices.second.x)
			<abs(edge.vertices.first.y-edge.vertices.second.y))
		{
			if(abs(edge.vertices.first.x-obc->GetGridBoundary(Left))
				<abs(edge.vertices.first.x-obc->GetGridBoundary(Right)))
				diff=Vector2D(obc->GetGridBoundary(Right)-obc->GetGridBoundary(Left),0);
			else
				diff=Vector2D(-obc->GetGridBoundary(Right)+obc->GetGridBoundary(Left),0);
		}
		else
		{
			if(abs(edge.vertices.first.y-obc->GetGridBoundary(Down))
				<abs(edge.vertices.first.y-obc->GetGridBoundary(Up)))
				diff=Vector2D(0,obc->GetGridBoundary(Up)-obc->GetGridBoundary(Down));
			else
				diff=Vector2D(0,-obc->GetGridBoundary(Up)+obc->GetGridBoundary(Down));
		}
		return diff;
	}

#ifdef RICH_MPI
	vector<int> GetCornerNeighbors(Tessellation const& v,int rank)
	{
		vector<int> neighbors=v.GetNeighbors(rank);
		sort(neighbors.begin(),neighbors.end());
		int n=static_cast<int>(neighbors.size());
		vector<int> minremove(1,-1);
		neighbors=RemoveList(neighbors,minremove);
		// Arrange the corners accordingly
		vector<int> edgeindex;
		ConvexEdges(edgeindex,&v,rank);
		vector<int> result(static_cast<size_t>(n));
		for(int i=0;i<n;++i)
		{
			int other=v.GetEdge(edgeindex[static_cast<size_t>(i)]).neighbors.first;
			if(other==rank)
				other=v.GetEdge(edgeindex[static_cast<size_t>(i)]).neighbors.second;
			int nextneigh=v.GetEdge(edgeindex[static_cast<size_t>((i+1)%n)]).neighbors.first;
			if(nextneigh==rank)
			  nextneigh=v.GetEdge(edgeindex[static_cast<size_t>((i+1)%n)]).neighbors.second;
			if(other==-1&&nextneigh==-1)
			{
				result[static_cast<size_t>(i)]=-1;
				continue;
			}
			if(other==-1||nextneigh==-1)
			{
				vector<int> nextedges;
				(nextneigh==-1)?ConvexEdges(nextedges,&v,other):ConvexEdges(nextedges,&v,nextneigh);
				int nedges=static_cast<int>(nextedges.size());
				int counter=0;
				for(int k=0;k<nedges;++k)
				{
					Edge e=v.GetEdge(nextedges[static_cast<size_t>(k)]);
					if(e.neighbors.first==rank||e.neighbors.second==rank)
					{
						counter=k;
						break;
					}
				}
				Edge nextedge=(other==-1)?v.GetEdge(nextedges[static_cast<size_t>((counter+2)%nedges)])
				  :v.GetEdge(nextedges[static_cast<size_t>((counter-2+nedges)%nedges)]);
				if(nextneigh==-1)
					result[static_cast<size_t>(i)]=(nextedge.neighbors.first==other)?
					nextedge.neighbors.second:nextedge.neighbors.first;
				else
					result[static_cast<size_t>(i)]=(nextedge.neighbors.first==nextneigh)?
					nextedge.neighbors.second:nextedge.neighbors.first;
			}
			else
			{
				// Do they have a mutual edge?
				vector<int> const& otheredges=v.GetCellEdges(other);
				int N=static_cast<int>(otheredges.size());
				int k;
				for(k=0;k<N;++k)
				{
					Edge const& edge=v.GetEdge(otheredges[static_cast<size_t>(k)]);
					if(edge.neighbors.first==nextneigh||
						edge.neighbors.second==nextneigh)
					{
						// We have a mutual edge
						Vector2D otherv = (abs(abs(v.GetMeshPoint(rank) - edge.vertices.first) - abs(v.GetMeshPoint(other) -
							edge.vertices.first))<abs(abs(v.GetMeshPoint(rank) - edge.vertices.second) - abs(v.GetMeshPoint(other) -
							edge.vertices.second))) ? edge.vertices.second : edge.vertices.first;
						// Find mutual neighbors
						vector<int> minremove2(1,-1);
						minremove2.push_back(rank);
						vector<int> neighbors1=v.GetNeighbors(other);
						vector<int> neighbors2=v.GetNeighbors(nextneigh);
						neighbors1=RemoveList(neighbors1,minremove2);
						neighbors2=RemoveList(neighbors2,minremove2);
						sort(neighbors1.begin(),neighbors1.end());
						sort(neighbors2.begin(),neighbors2.end());
						vector<int> restemp(10,-1);
						vector<int>::iterator it=set_intersection(neighbors1.begin(),
							neighbors1.end(),neighbors2.begin(),neighbors2.end(),restemp.begin());
						int size=static_cast<int>(it-restemp.begin());
						if(size==0)
						{
							result[static_cast<size_t>(i)]=-1;
							break;
						}
						if (size == 1)
						{
							const double d2 = abs(v.GetMeshPoint(other) - otherv);
							const double d1 = abs(abs(v.GetMeshPoint(restemp[0]) - otherv) - d2);
							const double eps2 = 1e-6;
							if (d1 > eps2*d2)
							{
								result[static_cast<size_t>(i)] = -1;
								break;
							}
						}
						vector<double> dist(static_cast<size_t>(size));
						for(int kk=0;kk<size;++kk)
							dist[static_cast<size_t>(kk)]=abs(otherv-v.GetMeshPoint(restemp[static_cast<size_t>(kk)]));
						result[static_cast<size_t>(i)]=restemp[static_cast<size_t>(min_element(dist.begin(),dist.end())
										   -dist.begin())];
						break;
					}
				}
				if(k<N)
					continue;
				// No mutual edge, probably a square, find mutual neighbor
				vector<int> neigh1=v.GetNeighbors(other);
				vector<int> neigh2=v.GetNeighbors(nextneigh);
				sort(neigh1.begin(),neigh1.end());
				sort(neigh2.begin(),neigh2.end());
				vector<int> minremove2(1,-1);
				minremove2.push_back(rank);
				neigh1=RemoveList(neigh1,minremove2);
				neigh2=RemoveList(neigh2,minremove2);
				vector<int> restemp(10,-1);
				set_intersection(neigh1.begin(),neigh1.end(),neigh2.begin(),
					neigh2.end(),restemp.begin());
				result[static_cast<size_t>(i)]=restemp[0];
			}
		}
		return result;
	}

			  /*
	int SumV(vector<vector<int> > const& v)
	{
		int res=0;
		int n=static_cast<int>(v.size());
		for(int i=0;i<n;++i)
		{
		  res+=static_cast<int>(v[static_cast<size_t>(i)].size());
		}
		return res;
	}
			  */
#endif // RICH_MPI

	template <typename T>
	bool EmptyVectorVector(vector<vector<T> > const& v)
	{
		int n=static_cast<int>(v.size());
		for(int i=0;i<n;++i)
		{
			if(!v[static_cast<size_t>(i)].empty())
				return false;
		}
		return true;
	}

	template <typename T>
	vector<vector<T> > CombineVectorVector(vector<vector<T> > const& v1,
		vector<vector<T> > const& v2)
	{
		vector<vector<T> > res(v1);
		assert(v1.size()==v2.size());
		int n=static_cast<int>(v1.size());
		for(int i=0;i<n;++i)
		{
			if(!v2[static_cast<size_t>(i)].empty())
				res[static_cast<size_t>(i)].insert(res[static_cast<size_t>(i)].end(),v2[static_cast<size_t>(i)].begin(),v2[static_cast<size_t>(i)].end());
		}
		return res;
	}
}

VoronoiMesh::VoronoiMesh(vector<Vector2D> const& points,OuterBoundary const& bc):
	logger(0),
	vertices_(),
	eps(1e-8),
	obc(0),
	cell_edges(vector<Edge> ()),
	edges(vector<Edge>()),
	CM(vector<Vector2D> ()),
	mesh_vertices(vector<vector<int> >()),
	Tri(),
	GhostProcs(vector<int> ()),
	GhostPoints(vector<vector<int> > ()),
	SentProcs(vector<int> ()),
	SentPoints(vector<vector<int> > ()),
	selfindex(vector<int> ()),
	NGhostReceived(vector<vector<int> > ()),
	OrgCorner(),
	Nextra(0)
{
	Initialise(points,&bc);
}

#ifdef RICH_MPI
VoronoiMesh::VoronoiMesh(vector<Vector2D> const& points,Tessellation const& proctess,
	OuterBoundary const& bc):
	logger(0),
	eps(1e-8),
	obc(0),
	cell_edges(vector<Edge> ()),
	edges(vector<Edge>()),
	CM(vector<Vector2D> ()),
	mesh_vertices(vector<vector<int> >()),
	Tri(),
	GhostProcs(vector<int> ()),
	GhostPoints(vector<vector<int> > ()),
	SentProcs(vector<int> ()),
	SentPoints(vector<vector<int> > ()),
	selfindex(vector<int> ()),
	NGhostReceived(vector<vector<int> > ()),
	OrgCorner(),
	Nextra(0)
{
	Initialise(points,proctess,&bc);
}
#endif

vector<int> VoronoiMesh::AddPointsAlongEdge(size_t point,vector<vector<int> > const&copied,
	int side)
{
  int ncopy=static_cast<int>(copied[static_cast<size_t>(side)].size());
	Vector2D vec=Tri.get_point(point);
	vector<double> dist(static_cast<size_t>(ncopy));
	for(size_t i=0;i<copied[static_cast<size_t>(side)].size();++i)
	  dist[i]=vec.distance(Tri.get_point(static_cast<size_t>(copied[static_cast<size_t>(side)][i])));
	const int copylength=min(7,static_cast<int>(copied[static_cast<size_t>(side)].size())-1);
	vector<int> index,toadd(static_cast<size_t>(copylength));
	sort_index(dist,index);
	for(int i=0;i<copylength;++i)
	  toadd[static_cast<size_t>(i)]=copied[static_cast<size_t>(side)][static_cast<size_t>(index[static_cast<size_t>(i)+1])];
	return toadd;
}

Vector2D VoronoiMesh::CalcFaceVelocity(Vector2D wl, Vector2D wr,Vector2D rL, Vector2D rR,
	Vector2D f)const
{
	const Vector2D wprime = ScalarProd(wl-wr,f-(rR+rL)/2)*(rR-rL)/pow(abs(rR-rL),2);
	return 0.5*(wl+wr) + wprime;
}

bool VoronoiMesh::NearBoundary(int index) const
{
	const int n=int(mesh_vertices[static_cast<size_t>(index)].size());
	const int N=Tri.get_length();
	for(int i=0;i<n;++i)
	{
	  const int n0=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].neighbors.first;
	  const int n1=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].neighbors.second;
		if(n0<0||n1<0||n0>=N||n1>=N)
			return true;
	}
	return false;
}

VoronoiMesh::~VoronoiMesh(void) {}

int VoronoiMesh::GetOriginalIndex(int point) const
{
	int npoints=GetPointNo();
	if(point<npoints)
		return point;
	else
	{
#ifdef RICH_MPI
		return point;
#else
		int counter=0;
		if(point<Nextra)
		{
			UniversalError eo("Tried to get original index of non exsistent cell");
			eo.AddEntry("Tried accessing cell",point);
			throw eo;
		}
		int maxcor=static_cast<int>(Tri.getCor().size());
		int cumulative=Nextra;
		while(cumulative<maxcor)
		{
		  int temp=static_cast<int>(GhostPoints[static_cast<size_t>(counter)].size());
			if((cumulative+temp)<=point)
			{
				cumulative+=temp;
				++counter;
			}
			else
			{
			  return GhostPoints[static_cast<size_t>(counter)][static_cast<size_t>(point-cumulative)];
			}
		}
		UniversalError eo("Tried to get original index of non exsistent cell");
		eo.AddEntry("Tried accessing cell",point);
		throw eo;
		#endif
	}
}

vector<int> VoronoiMesh::GetSelfPoint(void)const
{
	return selfindex;
}

VoronoiMesh::VoronoiMesh():
	logger(0),
	vertices_(),
	eps(1e-8),
	obc(0),
	cell_edges(vector<Edge> ()),
	edges(vector<Edge>()),
	CM(vector<Vector2D> ()),
	mesh_vertices(vector<vector<int> >()),
	Tri(),
	GhostProcs(vector<int> ()),
	GhostPoints(vector<vector<int> > ()),
	SentProcs(vector<int> ()),
	SentPoints(vector<vector<int> > ()),
	selfindex(vector<int> ()),
	NGhostReceived(vector<vector<int> > ()),
	OrgCorner(),
	Nextra(0) {}

VoronoiMesh::VoronoiMesh(VoronoiMesh const& other):
  logger(other.logger),
  vertices_(other.vertices_),
  eps(other.eps),
  obc(other.obc),
  cell_edges(other.cell_edges),
  edges(other.edges),
  CM(other.CM),
  mesh_vertices(other.mesh_vertices),
  Tri(other.Tri),
  GhostProcs(other.GhostProcs),
  GhostPoints(other.GhostPoints),
  SentProcs(other.SentProcs),
  SentPoints(other.SentPoints),
  selfindex(other.selfindex),
  NGhostReceived(other.NGhostReceived),
  OrgCorner(),
  Nextra(other.Nextra)
{}

void VoronoiMesh::build_v()
{
	Vector2D center,center_temp;
	int j;
	facet to_check;
	Edge edge_temp;
	Vector2D p_temp;
	mesh_vertices.clear();
	mesh_vertices.resize(static_cast<size_t>(Tri.get_length()));
	edges.reserve(static_cast<size_t>(Tri.get_length()*3.5));
	int N=Tri.GetOriginalLength();
	for(int i=0;i<N;++i)
		mesh_vertices[static_cast<size_t>(i)].reserve(7);
	int Nfacets=Tri.get_num_facet();
	vector<Vector2D> centers(static_cast<size_t>(Nfacets));
	for(int i=0;i<Nfacets;++i)
		centers[static_cast<size_t>(i)]=Tri.GetCircleCenter(i);
	for(int i=0;i<Nfacets;++i)
	{
		center=centers[static_cast<size_t>(i)];
		to_check=Tri.get_facet(i);
		for(j=0;j<3;++j)
		{
			if(to_check.neighbors[static_cast<size_t>(j)]==Tri.get_last_loc())
				continue;
			if(to_check.neighbors[static_cast<size_t>(j)]<i)
				continue;
			center_temp=centers[static_cast<size_t>(to_check.neighbors[static_cast<size_t>(j)])];
			{
				edge_temp.vertices.first = center;
				edge_temp.vertices.second = center_temp;
				edge_temp.neighbors.first = to_check.vertices[static_cast<size_t>(j)];
				edge_temp.neighbors.second = to_check.vertices[static_cast<size_t>(j+1)%3];

				if(legal_edge(&edge_temp))
				{
					// I added a change here, if edge has zero length I don't add it.
					if(edge_temp.GetLength()>eps*sqrt(Tri.GetFacetRadius(i)*
						Tri.GetFacetRadius(to_check.neighbors[static_cast<size_t>(j)])))
					{
						{
							if(edge_temp.neighbors.first<Tri.GetOriginalLength())
								mesh_vertices[static_cast<size_t>(edge_temp.neighbors.first)].push_back(static_cast<int>(edges.size()));
							if(edge_temp.neighbors.second<Tri.GetOriginalLength())
								mesh_vertices[static_cast<size_t>(edge_temp.neighbors.second)].push_back(static_cast<int>(edges.size()));
							edges.push_back(edge_temp);
						}
					}
				}
			}
		}
	}
}

namespace 
{
  vector<Vector2D> calc_procpoints(const OuterBoundary& bc)
  {
    vector<Vector2D> res(4);
    res[0] = Vector2D(bc.GetGridBoundary(Left),bc.GetGridBoundary(Down));
    res[1] = Vector2D(bc.GetGridBoundary(Right),bc.GetGridBoundary(Down));
    res[2] = Vector2D(bc.GetGridBoundary(Right),bc.GetGridBoundary(Up));
    res[3] = Vector2D(bc.GetGridBoundary(Left),bc.GetGridBoundary(Up));
    return res;
  }

  Vector2D GetReflection(OuterBoundary const& bc, size_t index, Vector2D const& point)
  {
	  switch (index)
	  {
	  case 0:
		  return point + Vector2D(2 * (bc.GetGridBoundary(Right) - point.x),0);
	  case 1:
		  return point + Vector2D(0 , 2 * (bc.GetGridBoundary(Up) - point.y));
	  case 2:
		  return point + Vector2D(2 * (bc.GetGridBoundary(Left) - point.x), 0);
	  case 3:
		  return point + Vector2D(0, 2 * (bc.GetGridBoundary(Down) - point.y));
	  }
	  throw UniversalError("Wrong index in VoronoiMesh::GetReflection");
  }
}

void VoronoiMesh::Initialise(vector<Vector2D>const& pv,OuterBoundary const* _bc)
{
	obc=_bc;
	Tri.build_delaunay(UpdatePoints(pv,obc),calc_procpoints(*obc));

	Nextra=static_cast<int>(Tri.ChangeCor().size());
	vector<vector<int> > toduplicate = Tri.BuildBoundary(_bc,_bc->GetBoxEdges());

	eps=1e-8;
	edges.clear();
	GhostPoints.clear();
	GhostProcs.clear();
	build_v();

	if(logger)
		logger->output(*this);

	CM.resize(Tri.getCor().size());
	for (size_t i = 0; i<pv.size(); ++i)
		CM[i] = CalcCellCM(i);

	size_t counter = pv.size() + 3;
	if(_bc->GetBoundaryType()==Periodic)
	{
		for(size_t i=0;i<8;++i)
		{
			GhostPoints.push_back(toduplicate[i]);
			GhostProcs.push_back(-1);
			for (size_t j = 0; j < toduplicate[i].size(); ++j)
			{
				CM[counter] = CM[toduplicate[i][j]] + (Tri.get_point(counter) -
					Tri.get_point(static_cast<size_t>(GetOriginalIndex(static_cast<int>(counter)))));
				++counter;
			}
		}
	}
	else
	{
		for (size_t i = 0; i < 4; ++i)
		{
			GhostPoints.push_back(toduplicate[i]);
			GhostProcs.push_back(-1);
			if (_bc->GetBoundaryType() == Rectengular()||(i%2)==1)
			{
				for (size_t j = 0; j < toduplicate[i].size(); ++j)
				{
					CM[counter] = GetReflection(*_bc, i, CM[toduplicate[i][j]]);
				}
			}
			else
			{
				for (size_t j = 0; j < toduplicate[i].size(); ++j)
				{
					CM[counter] = CM[toduplicate[i][j]] + (Tri.get_point(counter) -
						Tri.get_point(static_cast<size_t>(GetOriginalIndex(static_cast<int>(counter)))));
				}
			}
			++counter;
		}
	}
}

#ifdef RICH_MPI
void VoronoiMesh::Initialise(vector<Vector2D>const& pv,Tessellation const& vproc,
	OuterBoundary const* outer)
{
	NGhostReceived.clear();
	const int rank = get_mpi_rank();
	obc=outer;
	vector<int> cedges;
	ConvexEdges(cedges,&vproc,rank);
	cell_edges.clear();
	for(int i=0;i<static_cast<int>(cedges.size());++i)
		cell_edges.push_back(vproc.GetEdge(cedges[static_cast<size_t>(i)]));

	// Get the convex hull of the cell
	vector<Vector2D> cpoints;
	ConvexHull(cpoints,&vproc,rank);
	//Build the delaunay
	Tri.build_delaunay(pv,cpoints);
	eps=1e-8;
	edges.clear();
	GhostPoints.clear();
	GhostProcs.clear();
	NGhostReceived.clear();
	OrgCorner.clear();
	selfindex.resize(pv.size());
	size_t npv=pv.size();
	for(size_t i=0;i<npv;++i)
	  selfindex[i]=static_cast<int>(i);

	Nextra=static_cast<int>(Tri.ChangeCor().size());
	GhostPoints=Tri.BuildBoundary(outer,vproc,NGhostReceived,GhostProcs);
	build_v();

	vector<vector<int> > copied(GhostPoints),tempv;
	for(size_t i=0;i<copied.size();++i)
		if(!copied[i].empty())
			sort(copied[i].begin(),copied[i].end());
	GetCorners(copied,tempv);

	if(logger)
		logger->output(*this);
	
	CM.resize(pv.size());
	for(size_t i=0;i<pv.size();++i)
		CM[i]=CalcCellCM(i);

	// sort the points
	vector<int> indeces;
	sort_index(GhostProcs,indeces);
	sort(GhostProcs.begin(),GhostProcs.end());
	vector<vector<int> > temppoints,temppoints2;
	temppoints.push_back(GhostPoints[static_cast<size_t>(indeces[0])]);
	temppoints2.push_back(NGhostReceived[static_cast<size_t>(indeces[0])]);
	for(int i=1;i<static_cast<int>(GhostProcs.size());++i)
	{
		if(GhostProcs[static_cast<size_t>(i)]==GhostProcs[static_cast<size_t>(i)-1])
		{
			temppoints[temppoints.size()-1].insert(temppoints[temppoints.size()-1].end(),
				GhostPoints[static_cast<size_t>(indeces[static_cast<size_t>(i)])].begin(),GhostPoints[static_cast<size_t>(indeces[static_cast<size_t>(i)])].end());
			temppoints2[temppoints2.size()-1].insert(temppoints2[temppoints2.size()-1].end(),
				NGhostReceived[static_cast<size_t>(indeces[static_cast<size_t>(i)])].begin(),NGhostReceived[static_cast<size_t>(indeces[static_cast<size_t>(i)])].end());
		}
		else
		{
			temppoints.push_back(GhostPoints[static_cast<size_t>(indeces[static_cast<size_t>(i)])]);
			temppoints2.push_back(NGhostReceived[static_cast<size_t>(indeces[static_cast<size_t>(i)])]);
		}
		GhostProcs=unique(GhostProcs);
		NGhostReceived=temppoints2;
		GhostPoints=temppoints;
	}
	SendRecvVelocity(CM, GhostPoints,GhostProcs, NGhostReceived,Tri.getCor().size()));
}
#endif


bool VoronoiMesh::legal_edge(Edge *e) //checks if both ends of the edge are outside the grid and that the edge doesn't cross the grid
{
	if((e->neighbors.first<Tri.get_length())||
		(e->neighbors.second<Tri.get_length()))
		return true;
	else
		return false;
}

double VoronoiMesh::GetWidth(int index)const
{
	return sqrt(GetVolume(index)/M_PI);
}

vector<int> const& VoronoiMesh::GetCellEdges(int index) const
{
  return mesh_vertices.at(static_cast<size_t>(index));
}

double VoronoiMesh::GetVolume(int index) const
{
  const Vector2D center=Tri.get_point(static_cast<size_t>(index));
	double area=0;
	for (size_t i=0;i<mesh_vertices[static_cast<size_t>(index)].size();++i)
	{
	  const Vector2D p1 = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].vertices.first-
			center;
		const Vector2D p2 = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].vertices.second-
			center;
		area+=0.5*abs(ScalarProd(p1,zcross(p2)));
	}
	return area;
}

Vector2D VoronoiMesh::CalcCellCM(size_t index) const
{
  const Vector2D center=edges[static_cast<size_t>(mesh_vertices[index].front())].vertices.first;
	Vector2D pc(0,0);
	double area=0;
	for (size_t i=1;i<mesh_vertices[index].size();i++)
	{
	  const Edge& edge = edges[static_cast<size_t>(mesh_vertices[index][i])];
		const Vector2D p1 = edge.vertices.first - center;
		const Vector2D p2 = edge.vertices.second - center;
		const double area_temp = 0.5*abs(ScalarProd(p1,zcross(p2)));
		area += area_temp;
		pc += (area_temp/3.)*(center+edge.vertices.first+edge.vertices.second);
	}
	return pc/area;
}

vector<Vector2D>& VoronoiMesh::GetMeshPoints(void)
{
	return Tri.GetMeshPoints();
}

void VoronoiMesh::Update(vector<Vector2D> const& p)
{
	// Clean_up last step
	edges.clear();
	GhostPoints.clear();
	GhostProcs.clear();
	vector<Vector2D> procpoints;
	procpoints.push_back(Vector2D(obc->GetGridBoundary(Left),obc->GetGridBoundary(Down)));
	procpoints.push_back(Vector2D(obc->GetGridBoundary(Right),obc->GetGridBoundary(Down)));
	procpoints.push_back(Vector2D(obc->GetGridBoundary(Right),obc->GetGridBoundary(Up)));
	procpoints.push_back(Vector2D(obc->GetGridBoundary(Left),obc->GetGridBoundary(Up)));
	vector<Vector2D> points=UpdatePoints(p,obc);
	Tri.update(points,procpoints);

	Nextra=static_cast<int>(Tri.ChangeCor().size());
	vector<Edge> box_edges=obc->GetBoxEdges();
	vector<vector<int> > toduplicate=Tri.BuildBoundary(obc,box_edges);

	eps=1e-8;
	edges.clear();
	GhostPoints.clear();
	GhostProcs.clear();
	build_v();

	if(logger)
		logger->output(*this);

	CM.resize(Tri.getCor().size());
	for (size_t i = 0; i<p.size(); ++i)
		CM[i] = CalcCellCM(i);

	size_t counter = p.size() + 3;
	if(obc->GetBoundaryType()==Periodic)
	{
		for(size_t i=0;i<8;++i)
		{
			GhostPoints.push_back(toduplicate[i]);
			GhostProcs.push_back(-1);
			for (size_t j = 0; j < toduplicate[i].size(); ++j)
			{
				CM[counter] = CM[toduplicate[i][j]] + (Tri.get_point(counter) -
					Tri.get_point(static_cast<size_t>(GetOriginalIndex(static_cast<int>(counter)))));
				++counter;
			}
		}
	}
	else
	{
		for (size_t i = 0; i < 4; ++i)
		{
			GhostPoints.push_back(toduplicate[i]);
			GhostProcs.push_back(-1);
			if (obc->GetBoundaryType() == Rectengular()||(i%2)==1)
			{
				for (size_t j = 0; j < toduplicate[i].size(); ++j)
				{
					CM[counter] = GetReflection(*obc, i, CM[toduplicate[i][j]]);
				}
			}
			else
			{
				for (size_t j = 0; j < toduplicate[i].size(); ++j)
				{
					CM[counter] = CM[toduplicate[i][j]] + (Tri.get_point(counter) -
						Tri.get_point(static_cast<size_t>(GetOriginalIndex(static_cast<int>(counter)))));
				}
			}
			++counter;
		}
	}
}

#ifdef RICH_MPI
void VoronoiMesh::Update(vector<Vector2D> const& p,Tessellation const &vproc)
{
	// Clean_up last step
	edges.clear();
	GhostPoints.clear();
	GhostProcs.clear();
	selfindex.clear();
	NGhostReceived.clear();
	const int rank = get_mpi_rank();
	vector<int> cedges;
	ConvexEdges(cedges,&vproc,rank);
	cell_edges.clear();
	for(int i=0;i<static_cast<int>(cedges.size());++i)
		cell_edges.push_back(vproc.GetEdge(cedges[static_cast<size_t>(i)]));

	vector<Vector2D> cproc;
	ConvexHull(cproc,&vproc,rank);
	// Send/Recv points that entered and left domain
	vector<int> cornerproc=GetCornerNeighbors(vproc,rank);
	vector<int> proclist(vproc.GetCellEdges(rank).size());
	for(int i=0;i<static_cast<int>(proclist.size());++i)
		proclist[static_cast<size_t>(i)]=(cell_edges[static_cast<size_t>(i)].neighbors.first==rank)? vproc.GetOriginalIndex(
		cell_edges[static_cast<size_t>(i)].neighbors.second):vproc.GetOriginalIndex(cell_edges[static_cast<size_t>(i)].neighbors.first);
	vector<Vector2D> points=UpdateMPIPoints(vproc,rank,cornerproc,proclist,p,obc,
		selfindex,SentProcs,SentPoints);

	try
	{
		Tri.update(points,cproc);
	}
	catch(UniversalError &eo)
	{
		int NN=static_cast<int>(cproc.size());
		for(int j=0;j<NN;++j)
		{
			eo.AddEntry("cell point x",cproc[static_cast<size_t>(j)].x);
			eo.AddEntry("cell point y",cproc[static_cast<size_t>(j)].y);
		}
		eo.AddEntry("Error in rank",rank);
		voronoi_loggers::BinLogger log("verror"+int2str(rank)+".bin");
		log.output(vproc);
		throw;
	}
	Nextra=static_cast<int>(Tri.ChangeCor().size());
	GhostPoints=Tri.BuildBoundary(obc,vproc,NGhostReceived,GhostProcs);
	build_v();

	vector<vector<int> > copied(GhostPoints),tempv;
	for(size_t i=0;i<copied.size();++i)
		if(!copied[i].empty())
			sort(copied[i].begin(),copied[i].end());
	GetCorners(copied,tempv);


	if(logger)
		logger->output(*this);
	// sort the points
	vector<int> indeces;
	sort_index(GhostProcs,indeces);
	sort(GhostProcs.begin(),GhostProcs.end());
	vector<vector<int> > temppoints,temppoints2;
	temppoints.push_back(GhostPoints[static_cast<size_t>(indeces[0])]);
	temppoints2.push_back(NGhostReceived[static_cast<size_t>(indeces[0])]);
	for(int i=1;i<static_cast<int>(GhostProcs.size());++i)
	{
		if(GhostProcs[static_cast<size_t>(i)]==GhostProcs[static_cast<size_t>(i)-1])
		{
			temppoints[temppoints.size()-1].insert(temppoints[temppoints.size()-1].end(),
				GhostPoints[static_cast<size_t>(indeces[static_cast<size_t>(i)])].begin(),GhostPoints[static_cast<size_t>(indeces[static_cast<size_t>(i)])].end());
			temppoints2[temppoints2.size()-1].insert(temppoints2[temppoints2.size()-1].end(),
				NGhostReceived[static_cast<size_t>(indeces[static_cast<size_t>(i)])].begin(),NGhostReceived[static_cast<size_t>(indeces[static_cast<size_t>(i)])].end());
		}
		else
		{
			temppoints.push_back(GhostPoints[static_cast<size_t>(indeces[static_cast<size_t>(i)])]);
			temppoints2.push_back(NGhostReceived[static_cast<size_t>(indeces[static_cast<size_t>(i)])]);
		}
		GhostProcs=unique(GhostProcs);
		NGhostReceived=temppoints2;
		GhostPoints=temppoints;
	}
	CM.resize(p.size());
	for(size_t i=0;i<p.size();++i)
		CM[i] = CalcCellCM(i);
	SendRecvVelocity(CM, GhostPoints,GhostProcs, NGhostReceived,Tri.getCor().size()));
}
#endif

vector<int> VoronoiMesh::GetNeighbors(int index)const
{
  vector<int> res(mesh_vertices[static_cast<size_t>(index)].size());
  for(size_t i=0;i<res.size();++i)
    res[i] = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][i])].neighbors.first!=index ?
	edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][i])].neighbors.first :
	edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][i])].neighbors.second;
  return res;
}

vector<int> VoronoiMesh::GetLiteralNeighbors(int index)const
{
  int n=static_cast<int>(mesh_vertices[static_cast<size_t>(index)].size());
	vector<int> res;
	res.reserve(static_cast<size_t>(n));
	for(int i=0;i<n;++i)
	{
	  int other = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].neighbors.first;
		if(other!=index)
		{
			if(other>-1)
				res.push_back(other);
		}
		else
		{
			if(other>-1)
			  other=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].neighbors.second;
			res.push_back(other);
		}
	}
	return res;
}

Tessellation* VoronoiMesh::clone(void)const
{
	return new VoronoiMesh(*this);
}
/*
namespace
{
void CreateGetPrimitiveList(vector<int> const& ToRemove,vector<vector<int> >
const& Nghost,vector<vector<int> > &MPI_AMR_Send)
{
int nprocs=static_cast<int>(Nghost.size());
MPI_AMR_Send.resize(nprocs);
vector<vector<int> > SortedNghost(nprocs),SortIndex(nprocs);
// sort Nghost
for(int i=0;i<nprocs;++i)
{
SortedNghost[i]=Nghost[i];
sort_index(Nghost[i],SortIndex[i]);
sort(SortedNghost[i].begin(),SortedNghost[i].end());
}
for(int i=0;i<static_cast<int>(ToRemove.size());++i)
{
for(int j=0;j<nprocs;++j)
{
if(binary_search(SortedNghost[j].begin(),SortedNghost[j].end(),
ToRemove[i]))
{
int index2=lower_bound(SortedNghost[j].begin(),SortedNghost[j].end(),
ToRemove[i])-SortedNghost[j].begin();
MPI_AMR_Send[j].push_back(SortIndex[j][index2]);
}
}
}
}
}
*/

void VoronoiMesh::RemoveCells(vector<int> &ToRemove,vector<vector<int> > &VolIndex,
	vector<vector<double> > &Volratio)
{
	Remove_Cells(*this,ToRemove,VolIndex,Volratio);
}

void Remove_Cells(VoronoiMesh &V,vector<int> &ToRemove,
	vector<vector<int> > &VolIndex,vector<vector<double> > &Volratio)
{
#ifdef RICH_MPI
	vector<vector<int> > BoundaryRemove;
	vector<vector<vector<int> > > BoundaryNeigh;
	// Get rid of corner points
	vector<int> toremove2;
	for(int i=0;i<static_cast<int>(V.OrgCorner.size());++i)
	{
		if(V.OrgCorner[static_cast<size_t>(i)].empty())
			continue;
		for(int j=0;j<static_cast<int>(ToRemove.size());++j)
		{
			if(binary_search(V.OrgCorner[static_cast<size_t>(i)].begin(),V.OrgCorner[static_cast<size_t>(i)].end(),ToRemove[static_cast<size_t>(j)]))
				toremove2.push_back(ToRemove[static_cast<size_t>(j)]);
		}
	}
	sort(toremove2.begin(),toremove2.end());
	toremove2=unique(toremove2);
	ToRemove=RemoveList(ToRemove,toremove2);
	// BoundaryRemove is the list per proc what points are removed given by their indeces in the sent vector
	// BoundaryNeigh, for each point in boundary remove, what are the indeces in sentvector of the local neighbors
	V.FindBoundaryRemoveSend(ToRemove,BoundaryRemove,BoundaryNeigh);
	vector<vector<int> > LocalNeighbors,GhostNeighbors;
	SendRecvBoundaryRemove(BoundaryRemove,BoundaryNeigh,V,LocalNeighbors,
		GhostNeighbors);
	// First index in GhostNeighbors is the ghost point that is removed
	int Ntotal=V.GetPointNo();
#endif
	VolIndex.clear();
	Volratio.clear();
	size_t n=ToRemove.size();
	vector<int> RemovedEdges;
#ifdef RICH_MPI
	n+=LocalNeighbors.size();
#endif
	VolIndex.resize(n);
	Volratio.resize(n);
	for(size_t i=0;i<n;++i)
	{
		//double vol_inverse=1/V.GetVolume(ToRemove[i]);
		// Get neighboring points
		vector<Vector2D> neighcor;
		//Remove The boundary points
		vector<int> real_neigh; //Keeps track of the real points I copied
		real_neigh.reserve(10);
		vector<int> AllPoints;
#ifdef RICH_MPI
		if(i>=ToRemove.size())
		{
			real_neigh=LocalNeighbors[static_cast<size_t>(i)-ToRemove.size()];
			AllPoints=GhostNeighbors[static_cast<size_t>(i)-ToRemove.size()];
		}
#endif
		if(i<ToRemove.size())
		{
			vector<int> neigh=V.GetLiteralNeighbors(ToRemove[static_cast<size_t>(i)]);
			for(size_t j=0;j<neigh.size();++j)
			{
#ifdef RICH_MPI
				if(neigh[static_cast<size_t>(j)]>-1&&neigh[static_cast<size_t>(j)]<Ntotal)
					real_neigh.push_back(V.GetOriginalIndex(
					neigh[static_cast<size_t>(j)]));
				else
					if(neigh[static_cast<size_t>(j)]>=Ntotal)
						AllPoints.push_back(neigh[static_cast<size_t>(j)]);
#else
				if(neigh[static_cast<size_t>(j)]>-1)
					real_neigh.push_back(V.GetOriginalIndex(
					neigh[static_cast<size_t>(j)]));
#endif
			}
		}
		AllPoints.insert(AllPoints.end(),real_neigh.begin(),real_neigh.end());
#ifdef RICH_MPI
		vector<int> mpi_real_neigh=AllPoints;
		sort(mpi_real_neigh.begin(),mpi_real_neigh.end());
		mpi_real_neigh=unique(mpi_real_neigh);
#endif
		for(size_t j=0;j<real_neigh.size();++j)
		{
			if(real_neigh[static_cast<size_t>(j)]<V.GetPointNo())
			{
				vector<int> temp_neigh=V.GetNeighbors(real_neigh[static_cast<size_t>(j)]);
				for(size_t jj=0;jj<temp_neigh.size();++jj)
					if(temp_neigh[static_cast<size_t>(jj)]>-1)
						if(i>=ToRemove.size()||V.GetOriginalIndex(temp_neigh[static_cast<size_t>(jj)])!=ToRemove[static_cast<size_t>(i)])
							AllPoints.push_back(V.GetOriginalIndex(temp_neigh[static_cast<size_t>(jj)]));
			}
		}
		// save old volume
		vector<double> dv(real_neigh.size());
		vector<int> indeces;
		sort_index(real_neigh,indeces);
		sort(real_neigh.begin(),real_neigh.end());
		for(size_t j=0;j<real_neigh.size();++j)
			dv[static_cast<size_t>(j)]=V.GetVolume(real_neigh[static_cast<size_t>(j)]);
		sort(AllPoints.begin(),AllPoints.end());
		AllPoints=unique(AllPoints);
#ifdef RICH_MPI
		AllPoints=RemoveList(AllPoints,mpi_real_neigh);
#else
		AllPoints=RemoveList(AllPoints,real_neigh);
#endif
#ifdef RICH_MPI
		AllPoints.insert(AllPoints.begin(),mpi_real_neigh.begin(),mpi_real_neigh.end());
#else
		AllPoints.insert(AllPoints.begin(),real_neigh.begin(),real_neigh.end());
#endif
#ifdef RICH_MPI
		vector<int> remtemp;
		if(i>=ToRemove.size())
		{
			remtemp.push_back(GhostNeighbors[static_cast<size_t>(i)-ToRemove.size()][0]);
			AllPoints=RemoveList(AllPoints,remtemp);
			mpi_real_neigh=RemoveList(mpi_real_neigh,remtemp);
		}
#endif
		neighcor.reserve(AllPoints.size());
		for(size_t j=0;j<AllPoints.size();++j)
		  neighcor.push_back(V.Tri.get_point(static_cast<size_t>(AllPoints[static_cast<size_t>(j)])));
		// build new local tessellation
		VoronoiMesh Vlocal;
		Vlocal.Initialise(neighcor,V.obc);
		// save the volume change
		for(size_t j=0;j<real_neigh.size();++j)
		{
			double vtemp=Vlocal.GetVolume(static_cast<int>(j));
			dv[j]=(vtemp-dv[j]);
			if(dv[j]<(-vtemp*1e-6))
			{
				UniversalError eo("Negative volume difference in Remove cells");
				throw eo;
			}
		}
		Volratio[static_cast<size_t>(i)]=dv;
		VolIndex[static_cast<size_t>(i)]=real_neigh;
		// Fix the old Voronoi////////////////
		int Nedges=static_cast<int>(Vlocal.edges.size());
		// Keep only the relevant edges
		vector<Edge> NewEdges;
		int N=static_cast<int>(real_neigh.size());
		for(int j=0;j<Nedges;++j)
		{
			int temp0=Vlocal.edges[static_cast<size_t>(j)].neighbors.first;
			int temp1=Vlocal.edges[static_cast<size_t>(j)].neighbors.second;
			if(temp0<N)
			{
#ifndef RICH_MPI
				if(temp1<N||temp1<0||Vlocal.GetOriginalIndex(temp1)<N)
#else
				if(temp1<N||temp1<0||temp1<static_cast<int>(mpi_real_neigh.size()))
#endif
				{
					NewEdges.push_back(Vlocal.edges[static_cast<size_t>(j)]);
					continue;
				}
			}
			if(temp1<N)
#ifndef RICH_MPI
				if(temp0<0||Vlocal.GetOriginalIndex(temp0)<N)
#else
				if(temp0<0||temp0<static_cast<int>(mpi_real_neigh.size()))
#endif
					NewEdges.push_back(Vlocal.edges[static_cast<size_t>(j)]);
		}
		// Fix the new edges to correspond to the correct points in cor, use real_neigh
		for(size_t j=0;j<NewEdges.size();++j)
		{
			int temp=NewEdges[static_cast<size_t>(j)].neighbors.first;
			if(temp>=N)
			{
#ifndef RICH_MPI
				// New point add to tessellation
				V.Tri.AddAditionalPoint(Vlocal.GetMeshPoint(temp));
				NewEdges[static_cast<size_t>(j)].neighbors.first = V.Tri.GetCorSize()-1;
				// add it to the DuplicatedPoints
				V.GhostPoints[static_cast<size_t>(V.GhostPoints.size())-1].push_back(AllPoints[static_cast<size_t>
												  (Vlocal.GetOriginalIndex(temp))]);
#else
				NewEdges[static_cast<size_t>(j)].neighbors.first=AllPoints[static_cast<size_t>(NewEdges[static_cast<size_t>(j)].neighbors.first)];
#endif
			}
			else
			{
				if(temp>-1)
					NewEdges[static_cast<size_t>(j)].neighbors.first = real_neigh[static_cast<size_t>(temp)];
			}
			temp=NewEdges[static_cast<size_t>(j)].neighbors.second;
			if(temp>=N)
			{
#ifndef RICH_MPI
				// New point add to tessellation
				V.Tri.AddAditionalPoint(Vlocal.GetMeshPoint(temp));
				NewEdges[static_cast<size_t>(j)].neighbors.second = V.Tri.GetCorSize()-1;
				// add it to the DuplicatedPoints
				V.GhostPoints[static_cast<size_t>(V.GhostPoints.size())-1].push_back(AllPoints[static_cast<size_t>
												  (Vlocal.GetOriginalIndex(temp))]);
#else
				NewEdges[static_cast<size_t>(j)].neighbors.second=AllPoints[static_cast<size_t>(NewEdges[static_cast<size_t>(j)].neighbors.second)];
#endif
			}
			else
			{
				if(temp>-1)
					NewEdges[static_cast<size_t>(j)].neighbors.second = real_neigh[static_cast<size_t>(temp)];
			}
		}
		const int Npoints=V.GetPointNo();
		//Remove bad edge reference from mesh vertices and outer edges
		vector<int> oldedges;
		if(i<ToRemove.size())
		  oldedges=V.mesh_vertices[static_cast<size_t>(ToRemove[static_cast<size_t>(i)])];
		for(int j=0;j<N;++j)
		{
			int temp1=real_neigh[static_cast<size_t>(j)];
			int NN=static_cast<int>(V.mesh_vertices[static_cast<size_t>(temp1)].size());
			for(int jj=0;jj<NN;++jj)
			{
			  Edge etemp=V.edges[static_cast<size_t>(V.mesh_vertices[static_cast<size_t>(temp1)][static_cast<size_t>(jj)])];
#ifdef RICH_MPI
				if((i<ToRemove.size()&&(etemp.neighbors.first==-1||V.GetOriginalIndex(etemp.neighbors.first)
					==ToRemove[static_cast<size_t>(i)]))||((i>=ToRemove.size())&&(etemp.neighbors.first==-1||V.GetOriginalIndex(etemp.neighbors.first)
												      ==GhostNeighbors[static_cast<size_t>(i)-static_cast<size_t>(ToRemove.size())][0])))
#else
				if(etemp.neighbors.first==-1||V.GetOriginalIndex(etemp.neighbors.first)
					==ToRemove[i])
#endif
				{
					if(etemp.neighbors.first<0||etemp.neighbors.first>=Npoints)
						oldedges.push_back(V.mesh_vertices[static_cast<size_t>(temp1)][static_cast<size_t>(jj)]);
					V.mesh_vertices[static_cast<size_t>(temp1)].erase(V.mesh_vertices[static_cast<size_t>(temp1)].begin()
						+jj);
					--NN;
					--jj;
					continue;
				}
#ifdef RICH_MPI
				if(((i<ToRemove.size())&&(etemp.neighbors.second==-1||V.GetOriginalIndex(etemp.neighbors.second)
					==ToRemove[static_cast<size_t>(i)]))||((i>=ToRemove.size())&&(etemp.neighbors.second==-1||V.GetOriginalIndex(etemp.neighbors.second)
												      ==GhostNeighbors[static_cast<size_t>(i)-ToRemove.size()][0])))
#else
				if(etemp.neighbors.second==-1||V.GetOriginalIndex(etemp.neighbors.second)
					==ToRemove[static_cast<size_t>(i)])
#endif
				{
					if(etemp.neighbors.second<0||etemp.neighbors.second>=Npoints)
						oldedges.push_back(V.mesh_vertices[static_cast<size_t>(temp1)][static_cast<size_t>(jj)]);
					V.mesh_vertices[static_cast<size_t>(temp1)].erase(V.mesh_vertices[static_cast<size_t>(temp1)].begin()
						+jj);
					--NN;
					--jj;
					continue;
				}
			}
		}

		// copy the new edges
		Nedges=static_cast<int>(NewEdges.size());
		// Fix the mesh vertices to point to the correct edges
		int temp,other;
		for(int j=0;j<Nedges;++j)
		{
			const double R=NewEdges[static_cast<size_t>(j)].GetLength();
			const double eps=1e-8;
			for(int k=0;k<2;++k)
			{
				temp=pair_member(NewEdges[static_cast<size_t>(j)].neighbors,k);
				other=pair_member(NewEdges[static_cast<size_t>(j)].neighbors,(k+1)%2);
				if(temp==-1||temp>=Npoints)
					continue;
				size_t jj=0;
				if(other!=-1)
				{
					for(jj=0;jj<V.mesh_vertices[static_cast<size_t>(temp)].size();++jj)
					{
					  if(V.edges[static_cast<size_t>(V.mesh_vertices[static_cast<size_t>(temp)][static_cast<size_t>(jj)])].neighbors.first>-1)
						  if(V.GetMeshPoint(V.edges[static_cast<size_t>(V.mesh_vertices[static_cast<size_t>(temp)][static_cast<size_t>(jj)])].
								neighbors.first).distance(V.GetMeshPoint(other))<eps*R)
							{
							  V.edges[static_cast<size_t>(V.mesh_vertices[static_cast<size_t>(temp)][static_cast<size_t>(jj)])]=NewEdges[static_cast<size_t>(j)];
								break;
							}
						if(V.edges[static_cast<size_t>(V.mesh_vertices[static_cast<size_t>(temp)][static_cast<size_t>(jj)])].neighbors.second>-1)
							  if(V.GetMeshPoint(V.edges[static_cast<size_t>(V.mesh_vertices[static_cast<size_t>(temp)][static_cast<size_t>(jj)])].
									neighbors.second).distance(V.GetMeshPoint(other))<eps*R)
								{
								  V.edges[static_cast<size_t>(V.mesh_vertices[static_cast<size_t>(temp)][static_cast<size_t>(jj)])]=NewEdges[static_cast<size_t>(j)];
									break;
								}
					}
				}
				else
				{
					V.mesh_vertices[static_cast<size_t>(temp)].push_back(static_cast<int>(V.edges.size()));
					V.edges.push_back(NewEdges[static_cast<size_t>(j)]);
				}
				// We made New neighbors
				if(jj==V.mesh_vertices[static_cast<size_t>(temp)].size())
				{
					if(other<Npoints)
					{
					  V.mesh_vertices[static_cast<size_t>(temp)].push_back(static_cast<int>(V.edges.size())-k);
						if(k==0)
							V.edges.push_back(NewEdges[static_cast<size_t>(j)]);
					}
					else
					{
						V.mesh_vertices[static_cast<size_t>(temp)].push_back(static_cast<int>(V.edges.size()));
						V.edges.push_back(NewEdges[static_cast<size_t>(j)]);
					}
				}
			}
		}
		RemovedEdges.insert(RemovedEdges.end(),oldedges.begin(),oldedges.end());
	}

	// Done with all of the points to remove now start fixing
	sort(RemovedEdges.begin(),RemovedEdges.end());
	RemovedEdges=unique(RemovedEdges);
	sort(ToRemove.begin(),ToRemove.end());

	//Fix edges
	size_t Nedges=V.edges.size();
	int toReduce;
	for(size_t i=0;i<Nedges;++i)
	{
		int temp=V.edges[static_cast<size_t>(i)].neighbors.first;
		if(temp>-1)
		{
			toReduce=int(lower_bound(ToRemove.begin(),ToRemove.end(),temp)-
				ToRemove.begin());
			V.edges[static_cast<size_t>(i)].neighbors.first = temp - toReduce;
		}
		temp = V.edges[static_cast<size_t>(i)].neighbors.second;
		if(temp>-1)
		{
			toReduce=int(lower_bound(ToRemove.begin(),ToRemove.end(),temp)-
				ToRemove.begin());
			V.edges[static_cast<size_t>(i)].neighbors.second = temp - toReduce;
		}
	}
	// Remove bad edges
	RemoveVector(V.edges,RemovedEdges);
	// Fix mesh_vertices
	RemoveVector(V.mesh_vertices,ToRemove);
	size_t N=V.mesh_vertices.size();
	// Fix edge number in mesh_vertices
	for(size_t i=0;i<N;++i)
	{
		for(size_t j=0;j<V.mesh_vertices[static_cast<size_t>(i)].size();++j)
		{
			toReduce=int(lower_bound(RemovedEdges.begin(),RemovedEdges.end(),
				V.mesh_vertices[static_cast<size_t>(i)][static_cast<size_t>(j)])-RemovedEdges.begin());
			V.mesh_vertices[static_cast<size_t>(i)][static_cast<size_t>(j)]-=toReduce;
		}
	}

	// Fix cor in Tri
	RemoveVector(V.Tri.ChangeCor(),ToRemove);
	// Fix CM
	V.CM.resize(V.CM.size()-ToRemove.size());
	for(size_t i=0;i<V.CM.size();++i)
		V.CM[i]=V.CalcCellCM(i);
	n=ToRemove.size();
	// Fix point numer in Tri
	V.Tri.ChangeOlength(V.Tri.get_length()-static_cast<int>(n));
	V.Tri.Changelength(V.Tri.get_length()-static_cast<int>(n));
#ifdef RICH_MPI
	// Fix the self index
	for(int i=0;i<static_cast<int>(V.selfindex.size());++i)
	{
		int toReduce2=int(lower_bound(ToRemove.begin(),ToRemove.end(),V.selfindex[static_cast<size_t>(i)])-
			ToRemove.begin());
		V.selfindex[static_cast<size_t>(i)]-=toReduce2;
	}
	vector<int> ToRemovetemp;
	for(vector<int>::iterator it=ToRemove.begin();it!=ToRemove.end();++it)
		if(binary_search(V.selfindex.begin(),V.selfindex.end(),*it))
			ToRemovetemp.push_back(*it);
	RemoveList(V.selfindex,ToRemovetemp);
#endif
	// Fix the sent cells, do i need this??
	V.Nextra-=static_cast<int>(ToRemove.size());
	// Add the removed boundary with MPI
#ifdef RICH_MPI
	for(int i=0;i<static_cast<int>(GhostNeighbors.size());++i)
		if(!GhostNeighbors[static_cast<size_t>(i)].empty())
			ToRemove.push_back(GhostNeighbors[static_cast<size_t>(i)][0]);
#endif
	return;
	// Reset Tree if self gravity is needed
}

namespace
{
	int FindEdge(VoronoiMesh const& V,int tofind,int celltolook)
	{
		vector<int> edges=V.GetCellEdges(celltolook);
		int n=static_cast<int>(edges.size());
		for(int i=0;i<n;++i)
		{
			Edge edge=V.GetEdge(edges[static_cast<size_t>(i)]);
			if(V.GetOriginalIndex(edge.neighbors.first)==tofind||
				V.GetOriginalIndex(edge.neighbors.second)==tofind)
				return edges[static_cast<size_t>(i)];
		}
		throw("Couldn't find neighbor in Voronoi::FindEdge");
	}
}

int FixPeriodNeighbor(VoronoiMesh &V,int other,int ToRefine,int /*NewIndex*/,
	Vector2D const& NewPoint)
{
	int loc=FindEdge(V,ToRefine,other);
	vector<Vector2D>& cor=V.Tri.ChangeCor();
	cor.push_back(NewPoint);
	int& temp = V.edges[static_cast<size_t>(loc)].neighbors.second==other ?
		V.edges[static_cast<size_t>(loc)].neighbors.first :
	V.edges[static_cast<size_t>(loc)].neighbors.second;
	temp = static_cast<int>(cor.size());
	return static_cast<int>(cor.size());
}

void VoronoiMesh::RefineCells(vector<int> const& ToRefine,
	vector<Vector2D> const& directions,double alpha)
{
	Refine_Cells(*this,ToRefine,alpha,directions,obc->GetBoundaryType()==Periodic);
	return;
}

void Refine_Cells(VoronoiMesh &V,vector<int> const& ToRefine,double alpha,
		  vector<Vector2D> const& directions,
		  #ifndef RICH_MPI
		  bool PeriodicBoundary)
                  #else
                  bool /*PeriodicBoundary*/)
		  #endif
{
	int n=static_cast<int>(ToRefine.size());
	if(n==0)
		return;
	int Npoints=V.GetPointNo();
	vector<Vector2D>& cor=V.Tri.ChangeCor();
	vector<Vector2D> cortemp(cor.size()-static_cast<size_t>(Npoints));
	// copy the ghost points into cortemp
	copy(cor.begin()+Npoints,cor.end(),cortemp.begin());
	int N=static_cast<int>(cor.size());
	// Expand the cor vector to include the new points
	cor.resize(static_cast<size_t>(N+n));
	copy(cortemp.begin(),cortemp.end(),cor.begin()+Npoints+n);
	cortemp.clear();
	// Fix the boundary point refrences if needed

#ifndef RICH_MPI
	if(PeriodicBoundary)
	{
#endif
		for(int i=0;i<static_cast<int>(V.edges.size());++i)
		{
			int temp_neigh=V.edges[static_cast<size_t>(i)].neighbors.first;
			if(temp_neigh>Npoints)
			{
				V.edges[static_cast<size_t>(i)].neighbors.first = temp_neigh+n;
				continue;
			}
			temp_neigh=V.edges[static_cast<size_t>(i)].neighbors.second;
			if(temp_neigh>Npoints)
				V.edges[static_cast<size_t>(i)].neighbors.second = temp_neigh+n;
		}
#ifndef RICH_MPI
	}
#endif
	// Update the lengths
	V.Tri.Changelength(Npoints+n);
	V.Tri.ChangeOlength(Npoints+n);
	V.Nextra+=static_cast<int>(ToRefine.size());
	// reserve space for mesh_vertices
	V.mesh_vertices.resize(static_cast<size_t>(Npoints+n));
	// Refine the points
	for(int i=0;i<n;++i)
	{
		// Get the edges
		vector<int> edge_index=V.GetCellEdges(ToRefine[static_cast<size_t>(i)]);
		int nedges=static_cast<int>(edge_index.size());
		vector<Edge> edges(static_cast<size_t>(nedges));
		for(int j=0;j<nedges;++j)
			edges[static_cast<size_t>(j)]=V.GetEdge(edge_index[static_cast<size_t>(j)]);
		// Make the new point by finding the closest edge
		Vector2D NewPoint=V.GetMeshPoint(ToRefine[static_cast<size_t>(i)]);
		double R=V.GetWidth(ToRefine[static_cast<size_t>(i)]);
		Vector2D normal;
		Vector2D slope;
		// Do we have a specified direction or do we do line conneting the CM or nearest edge
		if(!directions.empty())
		{
			slope=directions[static_cast<size_t>(i)]/abs(directions[static_cast<size_t>(i)]);
			normal.Set(slope.y,-slope.x);
		}
		else
			slope=FindBestSplit(&V,ToRefine[static_cast<size_t>(i)],edges,R,normal);
		NewPoint+=alpha*R*slope;
		cor[static_cast<size_t>(Npoints+i)]=NewPoint;
		Vector2D v(V.GetMeshPoint(ToRefine[static_cast<size_t>(i)]));
		// Split edges and update neighbors
		vector<int> old_ref;
		vector<int> new_ref;
		// Calculate the splitting edge
		Edge splitedge(0.5*(NewPoint+V.GetMeshPoint(ToRefine[static_cast<size_t>(i)])) -
			normal*R*20,0.5*(NewPoint+V.GetMeshPoint(ToRefine[static_cast<size_t>(i)])) +
			normal*R*20,Npoints+i,ToRefine[static_cast<size_t>(i)]);
		Vector2D intersect;
		// Find other intersecting segment
		int counter_edges=0;
		int coincide;
		Vector2D vs1,vs2;
		for(int j=0;j<nedges;++j)
		{
			// Does the splitting edge cross the edge?
			if(SegmentIntersection(splitedge,edges[static_cast<size_t>(j)],intersect,1e-7))
			{
				if(counter_edges==0)
					vs1=intersect;
				else
					if(vs1.distance(intersect)<R*1e-8)
						continue;
					else
						vs2=intersect;
				++counter_edges;
				if(counter_edges==2)
				{
					splitedge.vertices.first=vs1;
					splitedge.vertices.second=vs2;
					break;
				}
			}
		}
		if(counter_edges!=2)
		{
			UniversalError eo("Error in splitting the cell");
			eo.AddEntry("Cell index",ToRefine[static_cast<size_t>(i)]);
			throw eo;
		}
		// add the splitting edge
		old_ref.push_back(static_cast<int>(V.edges.size()));
		new_ref.push_back(static_cast<int>(V.edges.size()));
		V.edges.push_back(splitedge);
		// make the new edges
		for(int j=0;j<nedges;++j)
		{
			// Does the edge have a coinciding point with the splitting edge?
			if(splitedge.vertices.first.distance(edges[static_cast<size_t>(j)].vertices.first)<
				1e-8*R||splitedge.vertices.first.distance(edges[static_cast<size_t>(j)].
				vertices.second)<1e-8*R)
				coincide=0;
			else
				if(splitedge.vertices.second.distance(edges[static_cast<size_t>(j)].vertices.first)<
					1e-8*R||splitedge.vertices.second.distance(edges[static_cast<size_t>(j)].
					vertices.second)<1e-8*R)
					coincide=1;
				else
					coincide=2;
			if(coincide<2)
			{
				// No need to change vertices only neighbors
				if(DistanceToEdge(v,edges[static_cast<size_t>(j)])>DistanceToEdge(NewPoint,edges[static_cast<size_t>(j)]))
				{
					Edge NewEdge(edges[static_cast<size_t>(j)]);
					int rindex=1;
					if(NewEdge.neighbors.first==ToRefine[static_cast<size_t>(i)])
					{
						rindex=0;
						NewEdge.neighbors.first=Npoints+i;
					}
					else
						NewEdge.neighbors.second=Npoints+i;
					V.edges[static_cast<size_t>(edge_index[static_cast<size_t>(j)])]=NewEdge;
					new_ref.push_back(edge_index[static_cast<size_t>(j)]);
					const Vector2D diff(V.GetMeshPoint(pair_member(NewEdge.neighbors,(rindex+1)%2))-V.GetMeshPoint(V.GetOriginalIndex(pair_member(NewEdge.neighbors,(rindex+1)%2))));
					if(pair_member(NewEdge.neighbors,(rindex+1)%2)>(n+Npoints))
						FixPeriodNeighbor(V,V.GetOriginalIndex(pair_member(NewEdge.neighbors,(rindex+1)%2)),
						ToRefine[static_cast<size_t>(i)],Npoints+i,NewPoint-diff);
				}
				else
					old_ref.push_back(edge_index[static_cast<size_t>(j)]);
				continue;
			}
			// Does the splitting edge cross the edge?
			SegmentIntersection(splitedge,edges[static_cast<size_t>(j)],intersect);
			if(DistanceToEdge(intersect,edges[static_cast<size_t>(j)])<1e-8*R)
			{
				Edge NewEdge;
				NewEdge.neighbors.first=Npoints+i;
				if(edges[static_cast<size_t>(j)].neighbors.first==ToRefine[static_cast<size_t>(i)])
					NewEdge.neighbors.second=edges[static_cast<size_t>(j)].neighbors.second;
				else
					NewEdge.neighbors.second=edges[static_cast<size_t>(j)].neighbors.first;
				int index;
				if(NewPoint.distance(edges[static_cast<size_t>(j)].vertices.first)<
					v.distance(edges[static_cast<size_t>(j)].vertices.first))
					index=0;
				else
					index=1;
				NewEdge.vertices.first = pair_member(edges[static_cast<size_t>(j)].vertices,index);
				NewEdge.vertices.second = intersect;
				Vector2D diff(0,0);
				if(NewEdge.neighbors.second>(n+Npoints))
				{
					diff=V.GetMeshPoint(NewEdge.neighbors.second)-
						V.GetMeshPoint(V.GetOriginalIndex(NewEdge.neighbors.second));
					Edge temp(NewEdge.vertices.first-diff,
						NewEdge.vertices.second-diff,
						Npoints+i,
						V.GetOriginalIndex(NewEdge.neighbors.second));
					V.mesh_vertices[static_cast<size_t>(temp.neighbors.second)].push_back(static_cast<int>(V.edges.size()));
					V.edges.push_back(temp);
					int loc=FindEdge(V,ToRefine[static_cast<size_t>(i)],temp.neighbors.second);
					int index2;
					if(temp.vertices.first.distance(V.edges[static_cast<size_t>(loc)].vertices.first)<
						temp.vertices.first.distance(V.edges[static_cast<size_t>(loc)].vertices.second))
						index2=0;
					else
						index2=1;
					set_pair_member(V.edges[static_cast<size_t>(loc)].vertices,index2,temp.vertices.second);
					FixPeriodNeighbor(V,V.GetOriginalIndex(NewEdge.neighbors.second),
						Npoints+i,Npoints+i,NewPoint-diff);
				}
				else
					if(NewEdge.neighbors.second>-1)
						V.mesh_vertices[static_cast<size_t>(NewEdge.neighbors.second)].push_back(static_cast<int>(V.edges.size()));
				new_ref.push_back(static_cast<int>(V.edges.size()));
				V.edges.push_back(NewEdge);
				// Do the other split
				NewEdge.vertices.first = pair_member(edges[static_cast<size_t>(j)].vertices,(index+1)%2);
				NewEdge.neighbors.first=ToRefine[static_cast<size_t>(i)];
				V.edges[static_cast<size_t>(edge_index[static_cast<size_t>(j)])]=NewEdge;
				old_ref.push_back(edge_index[static_cast<size_t>(j)]);
				continue;
			}
			// No need to split the edge
			// check if update is needed
			if(NewPoint.distance(edges[static_cast<size_t>(j)].vertices.first)<
				v.distance(edges[static_cast<size_t>(j)].vertices.first))
			{
				// Change edge neighbors
				int index=1;
				if(edges[static_cast<size_t>(j)].neighbors.first==ToRefine[static_cast<size_t>(i)])
				{
					index=0;
					V.edges[static_cast<size_t>(edge_index[static_cast<size_t>(j)])].neighbors.first=Npoints+i;
				}
				else
				  V.edges[static_cast<size_t>(edge_index[static_cast<size_t>(j)])].neighbors.second=Npoints+i;
				// add new reference
				new_ref.push_back(edge_index[static_cast<size_t>(j)]);
				if(pair_member(V.edges[static_cast<size_t>(edge_index[static_cast<size_t>(j)])].neighbors,(index+1)%2)>(n+Npoints))
				{
					int other=pair_member(edges[static_cast<size_t>(j)].neighbors,(index+1)%2);
					const Vector2D diff=V.GetMeshPoint(other)-V.GetMeshPoint(
						V.GetOriginalIndex(other));
					other=V.GetOriginalIndex(other);
					FixPeriodNeighbor(V,other,ToRefine[static_cast<size_t>(i)],Npoints+i,NewPoint-diff);
				}
			}
			else
				old_ref.push_back(edge_index[static_cast<size_t>(j)]);
		}
		V.mesh_vertices[static_cast<size_t>(Npoints+i)]=new_ref;
		V.mesh_vertices[static_cast<size_t>(ToRefine[static_cast<size_t>(i)])]=old_ref;
	}

	// Calculate the new CM
	for(int i=0;i<static_cast<int>(ToRefine.size());++i)
	  V.CM.push_back(V.CalcCellCM(static_cast<size_t>(Npoints+i)));
	for(int i=0;i<static_cast<int>(ToRefine.size());++i)
	  V.CM[static_cast<size_t>(ToRefine[static_cast<size_t>(i)])]=V.CalcCellCM(static_cast<size_t>(ToRefine[static_cast<size_t>(i)]));
	// Fix the self send
	for(int i=0;i<static_cast<int>(ToRefine.size());++i)
		V.selfindex.push_back(Npoints+i);
	// Fix the Nghost
	for(int i=0;i<static_cast<int>(V.NGhostReceived.size());++i)
	  for(int j=0;j<static_cast<int>(V.NGhostReceived[static_cast<size_t>(i)].size());++j)
			V.NGhostReceived[static_cast<size_t>(i)][static_cast<size_t>(j)]+=static_cast<int>(ToRefine.size());
	return;
	// Reset Tree if self gravity is needed
}

int VoronoiMesh::GetPointNo(void) const
{
	return Tri.get_length();
}

Vector2D VoronoiMesh::GetMeshPoint(int index) const
{
  return Tri.get_point(static_cast<size_t>(index));
}

int VoronoiMesh::GetTotalSidesNumber(void) const
{
	return static_cast<int>(edges.size());
}

const vector<Edge>& VoronoiMesh::getAllEdges(void) const
{
  return edges;
}

Edge const& VoronoiMesh::GetEdge(int index) const
{
	return edges[static_cast<size_t>(index)];
}

Vector2D const& VoronoiMesh::GetCellCM(int index) const
{
	return CM[static_cast<size_t>(index)];
}

vector<Vector2D>& VoronoiMesh::GetAllCM(void)
{
	return CM;
}

void VoronoiMesh::FindIntersectingOuterPoints(vector<Edge> const&box_edges,vector<vector<int> >
	&boxduplicate,vector<vector<int> > const&firstduplicated)
{
	int n=static_cast<int>(box_edges.size());
	boxduplicate.resize(static_cast<size_t>(n));
	int N=static_cast<int>(mesh_vertices.size());
	if(N<20)
	{
		for(int i=0;i<n;++i)
			for(int j=0;j<N;++j)
				boxduplicate[static_cast<size_t>(i)].push_back(j);
		return;
	}

	N=static_cast<int>(firstduplicated.size());
	for(int i=0;i<N;++i)
	{
	  n=static_cast<int>(firstduplicated[static_cast<size_t>(i)].size());
		for(int j=0;j<n;++j)
		{
			vector<int> temp=CellIntersectOuterBoundary(box_edges,firstduplicated[static_cast<size_t>(i)][static_cast<size_t>(j)]);
			int jj=static_cast<int>(temp.size());
			if(jj>0)
			{
				for(int k=0;k<jj;++k)
				  boxduplicate[static_cast<size_t>(temp[static_cast<size_t>(k)])].push_back(firstduplicated[static_cast<size_t>(i)][static_cast<size_t>(j)]);
			}
		}
	}
	n=static_cast<int>(box_edges.size());
	for(int i=0;i<n;++i)
	{
		if(!boxduplicate[static_cast<size_t>(i)].empty())
		{
			sort(boxduplicate[static_cast<size_t>(i)].begin(),boxduplicate[static_cast<size_t>(i)].end());
			boxduplicate[static_cast<size_t>(i)]=unique(boxduplicate[static_cast<size_t>(i)]);
		}
	}
}

void VoronoiMesh::FindIntersectingPoints(vector<Edge> const &box_edges,
	vector<vector<int> > &toduplicate)
{
	int n=static_cast<int>(box_edges.size());
	toduplicate.resize(static_cast<size_t>(n));
	int N=static_cast<int>(mesh_vertices.size());
	if(N<20)
	{
		for(int i=0;i<n;++i)
			for(int j=0;j<N;++j)
				toduplicate[static_cast<size_t>(i)].push_back(j);
		return;
	}

	for(int i=0;i<N;++i)
	{
		vector<int> temp;
		temp=CellIntersectBoundary(box_edges,i);
		int j=static_cast<int>(temp.size());
		if(j>0)
		{
			for(int k=0;k<j;++k)
			  toduplicate[static_cast<size_t>(temp[static_cast<size_t>(k)])].push_back(i);
		}
	}
	for(int i=0;i<n;++i)
	{
		if(!toduplicate[static_cast<size_t>(i)].empty())
		{
			sort(toduplicate[static_cast<size_t>(i)].begin(),toduplicate[static_cast<size_t>(i)].end());
			toduplicate[static_cast<size_t>(i)]=unique(toduplicate[static_cast<size_t>(i)]);
		}
	}
}

vector<int> VoronoiMesh::CellIntersectBoundary(vector<Edge> const&box_edges,int cell)
{
  int ncell=static_cast<int>(mesh_vertices[static_cast<size_t>(cell)].size());
	int nbox=static_cast<int>(box_edges.size());
	vector<int> res;
	Vector2D intersect;
	for(int i=0;i<ncell;++i)
	{
		for(int j=0;j<nbox;++j)
		{
		  if(SegmentIntersection(box_edges[static_cast<size_t>(j)],edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])],
				intersect))
				res.push_back(j);
		}
	}
	sort(res.begin(),res.end());
	res=unique(res);
	int nintersect=static_cast<int>(res.size());
	if(nintersect>1)
	{
		vector<Vector2D> cpoints;
		ConvexHull(cpoints,this,cell);
		for(int i=0;i<nbox;++i)
			if(PointInCell(cpoints,box_edges[static_cast<size_t>(i)].vertices.first)||
				PointInCell(cpoints,box_edges[static_cast<size_t>(i)].vertices.second))
				res.push_back(i);
		sort(res.begin(),res.end());
		res=unique(res);
	}
	return res;
}

vector<int> VoronoiMesh::CellIntersectOuterBoundary(vector<Edge> const&box_edges,int cell)
{
  int ncell=static_cast<int>(mesh_vertices[static_cast<size_t>(cell)].size());
	int nbox=static_cast<int>(box_edges.size());
	vector<int> res;
	Vector2D intersect;
	boost::array<Vector2D,3> tocheck;
	for(int i=0;i<ncell;++i)
	{
		for(int j=0;j<nbox;++j)
		{
		  if(SegmentIntersection(box_edges[static_cast<size_t>(j)],edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])],
				intersect))
			{
			  double r=sqrt(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].GetLength()*
					box_edges[static_cast<size_t>(j)].GetLength());
				double eps1=1e-7;
				if(abs(orient2d(TripleConstRef<Vector2D>
						(box_edges[static_cast<size_t>(j)].vertices.second-box_edges[static_cast<size_t>(j)].vertices.first,
						 edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.second-box_edges[static_cast<size_t>(j)].vertices.first,
						 edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.first-box_edges[static_cast<size_t>(j)].vertices.first)))
				   <r*r*eps1)
					continue;
				if(DistanceToEdge(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.first,
					box_edges[static_cast<size_t>(j)])<eps1*r)
					continue;
				if(DistanceToEdge(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.second,
					box_edges[static_cast<size_t>(j)])<eps1*r)
					continue;
				if((intersect.distance(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.first)
				    >eps1*r)&&(intersect.distance(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.second)
					>eps1*r))
					res.push_back(j);
			}
		}
	}
	sort(res.begin(),res.end());
	res=unique(res);
	return res;
}


bool VoronoiMesh::CloseToBorder(int point,int &border)
{
	int olength=Tri.GetOriginalLength();
	int n=static_cast<int>(mesh_vertices[static_cast<size_t>(point)].size());
	for(int i=0;i<n;++i)
	{
	  if(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.second==point)
		  border=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.first;
		else
		  border=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.second;
		if(border>olength)
			return true;
	}
	return false;
}

vector<int> VoronoiMesh::GetBorderingCells(vector<int> const& copied,
	vector<int> const& totest,int tocheck,vector<int> tempresult,int outer)
{
	int border,test;
	int olength=Tri.GetOriginalLength();
	tempresult.push_back(tocheck);
	sort(tempresult.begin(),tempresult.end());
	int n=static_cast<int>(mesh_vertices[static_cast<size_t>(tocheck)].size());
	for(int i=0;i<n;++i)
	{
	  if(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(tocheck)][static_cast<size_t>(i)])].neighbors.second==tocheck)
		  test=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(tocheck)][static_cast<size_t>(i)])].neighbors.first;
		else
		  test=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(tocheck)][static_cast<size_t>(i)])].neighbors.second;
		if(test>=olength)
			continue;
		if(test<0)
			continue;
		if(CloseToBorder(test,border))
			if(border==outer)
				if(!binary_search(copied.begin(),copied.end(),test)&&
					!binary_search(totest.begin(),totest.end(),test)&&
					!binary_search(tempresult.begin(),tempresult.end(),test))
					tempresult=GetBorderingCells(copied,totest,test,tempresult,outer);
	}
	return tempresult;
}

void VoronoiMesh::GetAdditionalBoundary(vector<vector<int> > &copied,
	vector<vector<int> > &neighbors,vector<vector<int> > &totest)
{
  int nsides=static_cast<int>(copied.size());
	// Get all the neighbors
	neighbors.clear();
	neighbors.resize(static_cast<size_t>(nsides));
	for(int i=0;i<nsides;++i)
	{
		sort(copied[static_cast<size_t>(i)].begin(),copied[static_cast<size_t>(i)].end());
		// look if there are boundary points neighbors
		int n=static_cast<int>(totest[static_cast<size_t>(i)].size());
		for(int j=0;j<n;++j)
		{
			if(totest[static_cast<size_t>(i)][static_cast<size_t>(j)]==-1)
				continue;
			vector<int> toadd;
			int outer=0;
			if(CloseToBorder(totest[static_cast<size_t>(i)][static_cast<size_t>(j)],outer))
				toadd=GetBorderingCells(copied[static_cast<size_t>(i)],totest[static_cast<size_t>(i)],totest[static_cast<size_t>(i)][static_cast<size_t>(j)],toadd,outer);
			int nn=static_cast<int>(toadd.size());
			for(int k=0;k<nn;++k)
				neighbors[static_cast<size_t>(i)].push_back(toadd[static_cast<size_t>(k)]);
		}
		sort(neighbors[static_cast<size_t>(i)].begin(),neighbors[static_cast<size_t>(i)].end());
		neighbors[static_cast<size_t>(i)]=unique(neighbors[static_cast<size_t>(i)]);
		neighbors[static_cast<size_t>(i)]=RemoveList(neighbors[static_cast<size_t>(i)],copied[static_cast<size_t>(i)]);
	}
}

void VoronoiMesh::GetRealNeighbor(vector<int> &result,int point) const
{
	result.reserve(7);
	int n=static_cast<int>(mesh_vertices[static_cast<size_t>(point)].size());
	int olength=Tri.GetOriginalLength();
	for(int i=0;i<n;++i)
	{
	  if(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.first==point)
		{
		  if(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.second>-1&&
			   edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.second<olength)
			  result.push_back(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.second);
		}
		else
		{
		  if(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.first>-1&&
			   edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.first<olength)
			  result.push_back(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.first);
		}
	}
	sort(result.begin(),result.end());
	result=unique(result);
}

void VoronoiMesh::GetNeighborNeighbors(vector<int> &result,int point) const
{
	vector<int> neigh=GetNeighbors(point);
	result.clear();
	result.reserve(25);
	//GetRealNeighbor(neigh,point);
	int n=static_cast<int>(neigh.size());
	for(int i=0;i<n;++i)
	{
	  if (neigh[static_cast<size_t>(i)] < 0)
			continue;
	  vector<int> temp=GetNeighbors(GetOriginalIndex(neigh[static_cast<size_t>(i)]));
		//GetRealNeighbor(temp,neigh[static_cast<size_t>(i)]);
		for(size_t j=0;j<temp.size();++j)
			result.push_back(GetOriginalIndex(temp[j]));
	}
	sort(result.begin(),result.end());
	result=unique(result);
	// Remove self point
	RemoveVal(result, point);
}

void VoronoiMesh::GetNeighborNeighborsMPI(vector<int> &result,int point)
{
	vector<int> neigh;
	result.clear();
	GetRealNeighbor(neigh,point);
	int n=static_cast<int>(neigh.size());
	for(int i=0;i<n;++i)
	{
		vector<int> temp;
		GetRealNeighbor(temp,neigh[static_cast<size_t>(i)]);
		int N=static_cast<int>(temp.size());
		for(int j=0;j<N;++j)
			result.push_back(temp[static_cast<size_t>(j)]);
	}
	sort(result.begin(),result.end());
	unique(result);
}

void VoronoiMesh::GetCorners(vector<vector<int> > &copied,
	vector<vector<int> > &result)
{
	// copied should be sorted already
	int nsides=static_cast<int>(copied.size());
	result.clear();
	OrgCorner.clear();
	OrgCorner.resize(static_cast<size_t>(nsides));
	result.resize(static_cast<size_t>(nsides));
	vector<vector<int> > toadd(static_cast<size_t>(nsides));
	for(int i=0;i<nsides;++i)
	{
	  int n=static_cast<int>(copied[static_cast<size_t>(i)].size());
		for(int j=0;j<n;++j)
		{
		  if(binary_search(copied[static_cast<size_t>((i+1)%nsides)].begin(),copied[static_cast<size_t>((i+1)%nsides)].end(),
				copied[static_cast<size_t>(i)][static_cast<size_t>(j)]))
			{
				vector<int> temp;
				GetNeighborNeighborsMPI(temp,copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
				result[static_cast<size_t>(i)].insert(result[static_cast<size_t>(i)].end(),temp.begin(),temp.end());
				temp=AddPointsAlongEdge(static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]),copied,i);
				toadd[static_cast<size_t>((i+1)%nsides)].insert(toadd[static_cast<size_t>((i+1)%nsides)].end(),temp.begin(),
					temp.end());
				temp=GetNeighbors(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
				for(vector<int>::iterator it=temp.begin();it!=temp.end();++it)
					if(*it<GetPointNo()&&*it>-1)
						OrgCorner[static_cast<size_t>(i)].push_back(*it);
				OrgCorner[static_cast<size_t>(i)].push_back(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
			}
		  if(binary_search(copied[static_cast<size_t>((i-1+nsides)%nsides)].begin(),copied[static_cast<size_t>((i-1+nsides)%nsides)].end(),
				copied[static_cast<size_t>(i)][static_cast<size_t>(j)]))
			{
				vector<int> temp;
				GetNeighborNeighborsMPI(temp,copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
				result[static_cast<size_t>((i-1+nsides)%nsides)].insert(result[static_cast<size_t>((i-1+nsides)%nsides)].end()
					,temp.begin(),temp.end());
				temp=AddPointsAlongEdge(static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]),copied,i);
				toadd[static_cast<size_t>((i-1+nsides)%nsides)].insert(toadd[static_cast<size_t>((i-1+nsides)%nsides)].end(),
					temp.begin(),temp.end());
				temp=GetNeighbors(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
				for(vector<int>::iterator it=temp.begin();it!=temp.end();++it)
					if(*it<GetPointNo()&&*it>-1)
					  OrgCorner[static_cast<size_t>((i-1+nsides)%nsides)].push_back(*it);
				OrgCorner[static_cast<size_t>((i-1+nsides)%nsides)].push_back(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
			}
		}
	}
	for(int i=0;i<nsides;++i)
	{
		copied[static_cast<size_t>(i)].insert(copied[static_cast<size_t>(i)].end(),toadd[static_cast<size_t>(i)].begin(),toadd[static_cast<size_t>(i)].end());
		sort(copied[static_cast<size_t>(i)].begin(),copied[static_cast<size_t>(i)].end());
		copied[static_cast<size_t>(i)]=unique(copied[static_cast<size_t>(i)]);
		sort(result[static_cast<size_t>(i)].begin(),result[static_cast<size_t>(i)].end());
		result[static_cast<size_t>(i)]=unique(result[static_cast<size_t>(i)]);
		if(!OrgCorner[static_cast<size_t>(i)].empty())
		{
			sort(OrgCorner[static_cast<size_t>(i)].begin(),OrgCorner[static_cast<size_t>(i)].end());
			OrgCorner[static_cast<size_t>(i)]=unique(OrgCorner[static_cast<size_t>(i)]);
		}
	}
}

void VoronoiMesh::GetToTest(vector<vector<int> > &copied,vector<vector<int> > &totest)
{
  int nsides=static_cast<int>(copied.size());
	int olength=Tri.GetOriginalLength();
	// sort the vectors
	for(int i=0;i<nsides;++i)
		sort(copied[static_cast<size_t>(i)].begin(),copied[static_cast<size_t>(i)].end());
	totest.resize(static_cast<size_t>(nsides));
	int test=0;
	for(int i=0;i<nsides;++i)
	{
		vector<int> totest2;
		int ncopy=static_cast<int>(copied[static_cast<size_t>(i)].size());
		for(int j=0;j<ncopy;++j)
		{
		  int n=static_cast<int>(mesh_vertices[static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)])].size());
			for(int k=0;k<n;++k)
			{
			  if(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)])][static_cast<size_t>(k)])].neighbors.first==
					copied[static_cast<size_t>(i)][static_cast<size_t>(j)])
				  test=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)])][static_cast<size_t>(k)])].neighbors.second;
				else
				  test=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)])][static_cast<size_t>(k)])].neighbors.first;
				if(test<olength)
					totest2.push_back(test);
			}
		}
		sort(totest2.begin(),totest2.end());
		totest2=unique(totest2);
		totest[static_cast<size_t>(i)]=totest2;
	}
}

vector<int> VoronoiMesh::FindEdgeStartConvex(int point)
{
  int n=static_cast<int>(mesh_vertices[static_cast<size_t>(point)].size());
	Vector2D min_point;
	int min_index=0,p_index;
	if(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(0)])].vertices.first.x<
	   edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][0])].vertices.second.x)
	{
	  min_point=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][0])].vertices.first;
		p_index=0;
	}
	else
	{
	  min_point=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][0])].vertices.second;
		p_index=1;
	}
	for(int i=1;i<n;++i)
	{
	  double R=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].GetLength();
		if(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.first.x<(min_point.x-R*eps))
		{
		  min_point=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.first;
			min_index=i;
			p_index=0;
		}
		if(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.second.x<(min_point.x-R*eps))
		{
		  min_point=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.second;
			min_index=i;
			p_index=1;
		}
		if(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.first.x<(min_point.x+R*eps)&&
		   edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.first.y<min_point.y)
		{
		  min_point=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.first;
			min_index=i;
			p_index=0;
		}

		if(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.second.x<(min_point.x+R*eps)&&
		   edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.second.y<min_point.y)
		{
		  min_point=edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.second;
			min_index=i;
			p_index=1;
		}
	}
	vector<int> res(2);
	res[0]=min_index;
	res[1]=p_index;
	return res;
}

void VoronoiMesh::ConvexEdgeOrder(void)
{
	int n=static_cast<int>(mesh_vertices.size());
	for(int i=0;i<n;++i)
	{
		double R=GetWidth(i);
		vector<int> min_index=FindEdgeStartConvex(i);
		int p_loc=min_index[1];
		int edge_loc=mesh_vertices[static_cast<size_t>(i)][static_cast<size_t>(min_index[0])];
		int nedges=static_cast<int>(mesh_vertices[static_cast<size_t>(i)].size());
		std::list<int> elist;
		for(int j=0;j<nedges;++j)
		{
			if(j!=min_index[0])
				elist.push_back(mesh_vertices[static_cast<size_t>(i)][static_cast<size_t>(j)]);
		}
		vector<int> new_order;
		new_order.reserve(static_cast<size_t>(nedges));
		new_order.push_back(edge_loc);
		for(int j=0;j<nedges;++j)
		{
			if(j==min_index[0])
				continue;
			int nlist=static_cast<int>(elist.size());
			std::list<int>::iterator it=elist.begin();
			for(int k=0;k<nlist;++k)
			{
			  double temp0=pair_member(edges[static_cast<size_t>(edge_loc)].vertices,(p_loc+1)%2).distance(edges[static_cast<size_t>(*it)].vertices.first);
				if(temp0<eps*R)
				{
					p_loc=0;
					edge_loc=*it;
					elist.erase(it);
					new_order.push_back(edge_loc);
					break;
				}
				double temp1=pair_member(edges[static_cast<size_t>(edge_loc)].vertices,(p_loc+1)%2).distance(
														edges[static_cast<size_t>(*it)].vertices.second);
				if(temp1<eps*R)
				{
					p_loc=1;
					edge_loc=*it;
					new_order.push_back(edge_loc);
					elist.erase(it);
					break;
				}
				++it;
			}
		}
		mesh_vertices[static_cast<size_t>(i)]=new_order;
	}
}

vector<Edge>& VoronoiMesh::GetAllEdges(void)
{
	return edges;
}

#ifdef RICH_MPI
void VoronoiMesh::NonSendBoundary(vector<int> &
	proclist,vector<vector<int> > & data,Tessellation const& v,
	vector<vector<int> > &totest,vector<Edge> const& boxedges)
{
	vector<int> newproclist;
	vector<vector<int> > newdata;
	vector<vector<int> > newtotest;
	int n=static_cast<int>(proclist.size());
	int nproc=v.GetPointNo();
	const int rank = get_mpi_rank();
	for(int i=0;i<n;++i)
	{
		if(proclist[static_cast<size_t>(i)]==-1)
		{
			RigidBoundaryPoints(data[static_cast<size_t>(i)],boxedges[static_cast<size_t>(i)]);
			continue;
		}
		if(proclist[static_cast<size_t>(i)]<nproc)
		{
			newproclist.push_back(proclist[static_cast<size_t>(i)]);
			newdata.push_back(data[static_cast<size_t>(i)]);
			newtotest.push_back(totest[static_cast<size_t>(i)]);
			continue;
		}
		if(v.GetOriginalIndex(proclist[static_cast<size_t>(i)])==rank)
		{
			PeriodicBoundaryPoints(data[static_cast<size_t>(i)],i);
			newproclist.push_back(proclist[static_cast<size_t>(i)]);
			newdata.push_back(data[static_cast<size_t>(i)]);
			newtotest.push_back(totest[static_cast<size_t>(i)]);
			continue;
		}
		// The other proc is periodic
		newproclist.push_back(proclist[static_cast<size_t>(i)]);
		newdata.push_back(data[static_cast<size_t>(i)]);
		newtotest.push_back(totest[static_cast<size_t>(i)]);
	}
	proclist=newproclist;
	data=newdata;
	totest=newtotest;
}
#endif

#ifdef MPI_VERSION
void VoronoiMesh::SendRecv(vector<int> const& procorder,vector<int> const&
	proclist,vector<vector<int> > &data)
{
	int n=static_cast<int>(procorder.size());
	int nlist=static_cast<int>(proclist.size());
	const int rank = get_mpi_rank();
	for(int i=0;i<n;++i)
	{
	  const int index=static_cast<int>(find(proclist.begin(),proclist.end(),procorder[static_cast<size_t>(i)])-proclist.begin());
		// Do we talk with this processor?
		if(index<nlist)
		{
			// Create send data
		  int nsend=static_cast<int>(data[static_cast<size_t>(index)].size());
		  vector<double> send(2*static_cast<size_t>(nsend));
			// Arrange the points to send in Hilbert order
			vector<Vector2D> cortemp=VectorValues(Tri.ChangeCor(),data[static_cast<size_t>(index)]);
			vector<int> order=HilbertOrder(cortemp,static_cast<int>(cortemp.size()));
			ReArrangeVector(data[static_cast<size_t>(index)],order);
			for(int j=0;j<nsend;++j)
			{
			  send[2*static_cast<size_t>(j)]=cortemp[static_cast<size_t>(order[static_cast<size_t>(j)])].x;
			  send[2*static_cast<size_t>(j)+1]=cortemp[static_cast<size_t>(order[static_cast<size_t>(j)])].y;
			}
			// Recv data
			MPI_Status status;
			vector<double> recv;
			int nrecv;
			if(rank<procorder[static_cast<size_t>(i)])
			{
				if(nsend>0)
					MPI_Send(&send[0],2*nsend,MPI_DOUBLE,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD);
				else
				{
					double temp=0;
					MPI_Send(&temp,1,MPI_DOUBLE,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD);
				}
				MPI_Probe(procorder[static_cast<size_t>(i)],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
				recv.resize(static_cast<size_t>(nrecv));
				int rtag=status.MPI_TAG;
				if(rtag==0)
					MPI_Recv(&recv[0],nrecv,MPI_DOUBLE,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD,&status);
				else
				{
					double temp=0;
					MPI_Recv(&temp,1,MPI_DOUBLE,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD,&status);
				}
			}
			else
			{
				MPI_Probe(procorder[static_cast<size_t>(i)],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
				recv.resize(static_cast<size_t>(nrecv));
				if(status.MPI_TAG==0)
					MPI_Recv(&recv[0],nrecv,MPI_DOUBLE,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD,&status);
				else
				{
					double temp=0;
					MPI_Recv(&temp,1,MPI_DOUBLE,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD,&status);
				}
				if(nsend>0)
					MPI_Send(&send[0],2*nsend,MPI_DOUBLE,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD);
				else
				{
					double temp=0;
					MPI_Send(&temp,1,MPI_DOUBLE,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD);
				}
			}
			nrecv/=2;
			vector<Vector2D> toadd(static_cast<size_t>(nrecv));
			for(int j=0;j<nrecv;++j)
			  toadd[static_cast<size_t>(j)]=Vector2D(recv[2*static_cast<size_t>(j)],recv[2*static_cast<size_t>(j)+1]);
			// Add the points
			if(!toadd.empty())
			{
				int Ntemp=Tri.GetTotalLength();
				for(int j=0;j<static_cast<int>(toadd.size());++j)
					NGhostReceived[static_cast<size_t>(index)].push_back(Ntemp+j);
				Tri.AddBoundaryPoints(toadd);
			}
		}
	}
}
#endif

#ifdef MPI_VERSION
void VoronoiMesh::SendRecvRemove(vector<int> const& procorder,vector<int> const&
	proclist,vector<vector<int> > &data)
{
	int n=static_cast<int>(procorder.size());
	int nlist=static_cast<int>(proclist.size());
	const int rank = get_mpi_rank();
	for(int i=0;i<n;++i)
	{
	  const int index=static_cast<int>(find(proclist.begin(),proclist.end(),procorder[static_cast<size_t>(i)])-proclist.begin());
		// Do we talk with this processor?
		if(index<nlist)
		{
			GhostPoints[static_cast<size_t>(index)].insert(GhostPoints[static_cast<size_t>(index)].end(),data[static_cast<size_t>(index)].begin(),
				data[static_cast<size_t>(index)].end());
			// Create send data
			int nsend=static_cast<int>(data[static_cast<size_t>(index)].size());
			vector<double> send(2*static_cast<size_t>(nsend));
			vector<Vector2D> cortemp=VectorValues(Tri.ChangeCor(),data[static_cast<size_t>(index)]);
			for(int j=0;j<nsend;++j)
			{
			  send[2*static_cast<size_t>(j)]=cortemp[static_cast<size_t>(j)].x;
			  send[2*static_cast<size_t>(j)+1]=cortemp[static_cast<size_t>(j)].y;
			}
			// Recv data
			MPI_Status status;
			vector<double> recv;
			int nrecv;
			if(rank<procorder[static_cast<size_t>(i)])
			{
				if(nsend>0)
					MPI_Send(&send[0],2*nsend,MPI_DOUBLE,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD);
				else
				{
					double temp=0;
					MPI_Send(&temp,1,MPI_DOUBLE,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD);
				}
				MPI_Probe(procorder[static_cast<size_t>(i)],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
				recv.resize(static_cast<size_t>(nrecv));
				int rtag=status.MPI_TAG;
				if(rtag==0)
					MPI_Recv(&recv[0],nrecv,MPI_DOUBLE,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD,&status);
				else
				{
					double temp=0;
					MPI_Recv(&temp,1,MPI_DOUBLE,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD,&status);
				}
			}
			else
			{
				MPI_Probe(procorder[static_cast<size_t>(i)],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
				recv.resize(static_cast<size_t>(nrecv));
				if(status.MPI_TAG==0)
					MPI_Recv(&recv[0],nrecv,MPI_DOUBLE,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD,&status);
				else
				{
					double temp=0;
					MPI_Recv(&temp,1,MPI_DOUBLE,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD,&status);
				}
				if(nsend>0)
					MPI_Send(&send[0],2*nsend,MPI_DOUBLE,procorder[static_cast<size_t>(i)],0,MPI_COMM_WORLD);
				else
				{
					double temp=0;
					MPI_Send(&temp,1,MPI_DOUBLE,procorder[static_cast<size_t>(i)],1,MPI_COMM_WORLD);
				}
			}
			nrecv/=2;
			vector<Vector2D> toadd(static_cast<size_t>(nrecv));
			for(int j=0;j<nrecv;++j)
			  toadd[static_cast<size_t>(j)]=Vector2D(recv[2*static_cast<size_t>(j)],recv[2*static_cast<size_t>(j)+1]);
			// Add the points
			if(!toadd.empty())
			{
				int Ntemp=Tri.GetTotalLength();
				for(int j=0;j<static_cast<int>(toadd.size());++j)
					NGhostReceived[static_cast<size_t>(index)].push_back(Ntemp+j);
				vector<Vector2D>& cor=Tri.ChangeCor();
				cor.insert(cor.end(),toadd.begin(),toadd.end());
			}
		}
	}
}
#endif

void VoronoiMesh::RigidBoundaryPoints(vector<int> &points,Edge const& edge)
{
	int npoints=static_cast<int>(points.size());
	vector<Vector2D> toadd;
	vector<int> pointstemp;
	pointstemp.reserve(static_cast<size_t>(npoints));
	toadd.reserve(static_cast<size_t>(npoints));
	Vector2D par(Parallel(edge));
	par=par/abs(par);
	Vector2D edge0=edge.vertices.first;
	boost::array<double,4> maxedges=FindMaxCellEdges();
	double dx=maxedges[1]-maxedges[0];
	double dy=maxedges[3]-maxedges[static_cast<size_t>(2)];
	for(int i=0;i<npoints;++i)
	{
	  Vector2D point=Tri.get_point(static_cast<size_t>(points[static_cast<size_t>(i)]));
		Vector2D temp=point-edge0;
		temp=2*par*ScalarProd(par,temp)-temp+edge0;
		if((abs(point.x-temp.x)<2*dx)&&(abs(point.y-temp.y)<2*dy))
		{
			toadd.push_back(temp);
			pointstemp.push_back(i);
		}
	}
	if(!toadd.empty())
	{
		vector<int> order=HilbertOrder(toadd,static_cast<int>(toadd.size()));
		ReArrangeVector(toadd,order);
		Tri.AddBoundaryPoints(toadd);
		ReArrangeVector(pointstemp,order);
		points=pointstemp;
	}
}

void VoronoiMesh::PeriodicBoundaryPoints(vector<int> &points,int edge_number)
{
	int npoints=static_cast<int>(points.size());
	vector<Vector2D> toadd(static_cast<size_t>(npoints));
	Vector2D diff=GetPeriodicDiff(cell_edges[static_cast<size_t>(edge_number)],obc);
	for(int i=0;i<npoints;++i)
	  toadd[static_cast<size_t>(i)]=Tri.get_point(static_cast<size_t>(points[static_cast<size_t>(i)]))+diff;
	if(!toadd.empty())
	{
		vector<int> order=HilbertOrder(toadd,static_cast<int>(toadd.size()));
		ReArrangeVector(toadd,order);
		Tri.AddBoundaryPoints(toadd);
		ReArrangeVector(points,order);
	}
}

#ifdef RICH_MPI
void VoronoiMesh::NonSendCorners(vector<int> &
	proclist,vector<vector<int> > & data,Tessellation const& v)
{
	vector<int> newproclist;
	vector<vector<int> > newdata;
	int n=static_cast<int>(proclist.size());
	int nproc=v.GetPointNo();
	const int rank = get_mpi_rank();
	for(int i=0;i<n;++i)
	{
		if(proclist[static_cast<size_t>(i)]<nproc)
		{
			newproclist.push_back(proclist[static_cast<size_t>(i)]);
			newdata.push_back(data[static_cast<size_t>(i)]);
			continue;
		}
		if(v.GetOriginalIndex(proclist[static_cast<size_t>(i)])==rank)
		{
			CornerBoundaryPoints(data[static_cast<size_t>(i)],i);
			continue;
		}
		// The other proc is periodic
		newproclist.push_back(proclist[static_cast<size_t>(i)]);
		newdata.push_back(data[static_cast<size_t>(i)]);
	}
	proclist=newproclist;
	data=newdata;
}
#endif

void VoronoiMesh::CornerBoundaryPoints(vector<int> &points,int edge_number)
{
	int n=static_cast<int>(cell_edges.size());
	Vector2D diff1=GetPeriodicDiff(cell_edges[static_cast<size_t>(edge_number)],obc);
	Vector2D diff2=GetPeriodicDiff(cell_edges[static_cast<size_t>(((edge_number+1)%n))],obc);
	int npoints=static_cast<int>(points.size());
	vector<Vector2D> toadd(static_cast<size_t>(npoints));
	for(int i=0;i<npoints;++i)
	  toadd[static_cast<size_t>(i)]=Tri.get_point(static_cast<size_t>(points[static_cast<size_t>(i)]))+diff1+diff2;
	if(!toadd.empty())
	{
		vector<int> order=HilbertOrder(toadd,static_cast<int>(toadd.size()));
		ReArrangeVector(toadd,order);
		Tri.AddBoundaryPoints(toadd);
		ReArrangeVector(points,order);
	}
}

vector<vector<int> >& VoronoiMesh::GetDuplicatedPoints(void)
{
	return GhostPoints;
}

vector<vector<int> >const& VoronoiMesh::GetDuplicatedPoints(void)const
{
	return GhostPoints;
}

vector<vector<int> >& VoronoiMesh::GetGhostIndeces(void)
{
	return NGhostReceived;
}

vector<vector<int> >const& VoronoiMesh::GetGhostIndeces(void)const
{
	return NGhostReceived;
}

int VoronoiMesh::GetTotalPointNumber(void)const
{
  return static_cast<int>(Tri.getCor().size());
}

vector<int> VoronoiMesh::GetDuplicatedProcs(void)const
{
	return GhostProcs;
}

vector<int> VoronoiMesh::GetSentProcs(void)const
{
	return SentProcs;
}

vector<vector<int> >const& VoronoiMesh::GetSentPoints(void)const
{
	return SentPoints;
}

// cpoints must be convex hull, checks if vec is inside cpoints
bool PointInCell(vector<Vector2D> const& cpoints,Vector2D const& vec)
{
	for(size_t i=0, endp=cpoints.size();i<endp;++i)
	{
	  if(orient2d(TripleConstRef<Vector2D>(cpoints[i],
					       cpoints[(i+1)%endp],
					       vec))<0)
	  return false;
	}
	return true;
}

// result is : minx, maxx, miny, maxy
boost::array<double,4> VoronoiMesh::FindMaxCellEdges(void)
{
	int n=static_cast<int>(cell_edges.size());
	boost::array<double,4> res;
	res[0]=min(cell_edges[0].vertices.first.x,cell_edges[0].vertices.second.x);
	res[1]=max(cell_edges[0].vertices.first.x,cell_edges[0].vertices.second.x);
	res[2]=min(cell_edges[0].vertices.first.y,cell_edges[0].vertices.second.y);
	res[3]=max(cell_edges[0].vertices.first.y,cell_edges[0].vertices.second.y);
	for(int i=1;i<n;++i)
	{
		res[0]=min(min(cell_edges[static_cast<size_t>(i)].vertices.first.x,cell_edges[static_cast<size_t>(i)].vertices.second.x),res[0]);
		res[1]=max(max(cell_edges[static_cast<size_t>(i)].vertices.first.x,cell_edges[static_cast<size_t>(i)].vertices.second.x),res[1]);
		res[2]=min(min(cell_edges[static_cast<size_t>(i)].vertices.first.y,cell_edges[static_cast<size_t>(i)].vertices.second.y),res[2]);
		res[3]=max(max(cell_edges[static_cast<size_t>(i)].vertices.first.y,cell_edges[static_cast<size_t>(i)].vertices.second.y),res[3]);
	}
	return res;
}

// BoundaryRemove is the list per proc what points are removed given by their indeces in the sent vector
// BoundaryNeigh, for each point in boundary remove, what are the indeces in sentvector of the local neighbors
#ifdef RICH_MPI
void VoronoiMesh::FindBoundaryRemoveSend(vector<int> const& ToRemove,
	vector<vector<int> > &BoundaryRemove,
	vector<vector<vector<int> > > &BoundaryNeigh)
{
	const int rank = get_mpi_rank();
	const int ws = get_mpi_size();
	vector<int> procorder=GetProcOrder(rank,ws);

	int npoints=GetPointNo();
	int nremove=static_cast<int>(ToRemove.size());
	int nprocs=static_cast<int>(GhostProcs.size());
	vector<vector<int> > NewSend;
	BoundaryRemove.clear();
	NewSend.clear();
	BoundaryNeigh.clear();
	BoundaryRemove.resize(static_cast<size_t>(nprocs));
	NewSend.resize(static_cast<size_t>(nprocs));
	BoundaryNeigh.resize(static_cast<size_t>(nprocs));
	vector<vector<int> > sort_indeces(static_cast<size_t>(nprocs)),sort_indecesg(static_cast<size_t>(nprocs));
	vector<vector<int> > SortedNGhostPoints(static_cast<size_t>(nprocs)),SortedGhostPoints(static_cast<size_t>(nprocs));
	for(int i=0;i<nprocs;++i)
	{
		sort_index(NGhostReceived[static_cast<size_t>(i)],sort_indeces[static_cast<size_t>(i)]);
		SortedNGhostPoints[static_cast<size_t>(i)]=NGhostReceived[static_cast<size_t>(i)];
		sort(SortedNGhostPoints[static_cast<size_t>(i)].begin(),SortedNGhostPoints[static_cast<size_t>(i)].end());
		sort_index(GhostPoints[static_cast<size_t>(i)],sort_indecesg[static_cast<size_t>(i)]);
		SortedGhostPoints[static_cast<size_t>(i)]=GhostPoints[static_cast<size_t>(i)];
		sort(SortedGhostPoints[static_cast<size_t>(i)].begin(),SortedGhostPoints[static_cast<size_t>(i)].end());
	}
	for(int i=0;i<nremove;++i)
	{
		// Is the point near a boundary?
		vector<int> neigh=GetNeighbors(ToRemove[static_cast<size_t>(i)]);
		for(vector<int>::iterator it=neigh.begin();it!=neigh.end();++it)
		{
			if(*it>=npoints)
			{
				// We have a boundary point, get its local neighbors;
				vector<int> localneigh,ghostneigh;
				for(vector<int>::iterator it2=neigh.begin();it2!=neigh.end();++it2)
					if(*it2<npoints&&*it2>-1)
						localneigh.push_back(*it2);
					else
						if(*it>=npoints)
							ghostneigh.push_back(*it2);
				// To which cpu did we send this point?
				for(int j=0;j<nprocs;++j)
				{
					for(size_t ii=0;ii<ghostneigh.size();++ii)
					{
						if(binary_search(SortedNGhostPoints[static_cast<size_t>(j)].begin(),SortedNGhostPoints[static_cast<size_t>(j)].end(),
							ghostneigh[static_cast<size_t>(ii)]))
						{
						  int index=static_cast<int>(lower_bound(SortedGhostPoints[static_cast<size_t>(j)].begin(),
								SortedGhostPoints[static_cast<size_t>(j)].end(),ToRemove[static_cast<size_t>(i)])
									     -SortedGhostPoints[static_cast<size_t>(j)].begin());
							BoundaryRemove[static_cast<size_t>(j)].push_back(sort_indecesg[static_cast<size_t>(j)][static_cast<size_t>(index)]);
							// What are the indeces of the sent points?
							// Is there a local point not already sent?
							int nlocal=static_cast<int>(localneigh.size());
							vector<int> temp;
							for(int k=0;k<nlocal;++k)
							{
								if(!binary_search(SortedGhostPoints[static_cast<size_t>(j)].begin(),
									SortedGhostPoints[static_cast<size_t>(j)].end(),localneigh[static_cast<size_t>(k)]))
								{
									if(NewSend[static_cast<size_t>(j)].empty())
									{
										NewSend[static_cast<size_t>(j)].push_back(localneigh[static_cast<size_t>(k)]);
										temp.push_back(static_cast<int>(SortedGhostPoints[static_cast<size_t>(j)].size()));
									}
									else
									{
									  int index2=static_cast<int>(find(NewSend[static_cast<size_t>(j)].begin(),NewSend[static_cast<size_t>(j)].end(),
													  localneigh[static_cast<size_t>(k)])-NewSend[static_cast<size_t>(j)].begin());
										if(index2>=static_cast<int>(NewSend[static_cast<size_t>(j)].size()))
											NewSend[static_cast<size_t>(j)].push_back(localneigh[static_cast<size_t>(k)]);
										temp.push_back(static_cast<int>(SortedGhostPoints[static_cast<size_t>(j)].size())+index2);
									}
								}
								else
								{
								  int index2=static_cast<int>(lower_bound(SortedGhostPoints[static_cast<size_t>(j)].begin(),
													 SortedGhostPoints[static_cast<size_t>(j)].end(),localneigh[static_cast<size_t>(k)])-SortedGhostPoints[static_cast<size_t>(j)].begin());
									temp.push_back(sort_indecesg[static_cast<size_t>(j)][static_cast<size_t>(index2)]);
								}
							}
							BoundaryNeigh[static_cast<size_t>(j)].push_back(temp);
							break;
						}
					}
				}
				break;
			}
		}
	}
	SendRecvRemove(procorder,GhostProcs,NewSend);
}
#endif
