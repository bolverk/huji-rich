#include "RefineStrategy.hpp"

RefineStrategy::RefineStrategy(void):
refined_old(vector<int>()) {}

vector<int> RefineStrategy::RemoveNearBoundary(vector<int> const& ToRefine,vector<Vector2D>
	&directions,Tessellation
	const& tess)
{
	int nrefine=(int)ToRefine.size();
	int npoints=tess.GetPointNo();
	vector<int> res;
	vector<Vector2D> newdirections;
	for(int i=0;i<nrefine;++i)
	{
		vector<int> const& edges=tess.GetCellEdges(ToRefine[i]);
		int nedge=(int) edges.size();
		bool good=true;
		for(int j=0;j<nedge;++j)
		{
			Edge const& edge=tess.GetEdge(edges[j]);
			if(edge.GetNeighbor(0)>npoints||edge.GetNeighbor(1)>npoints)
			{
				good=false;
				break;
			}
		}
		if(good)
		{
			res.push_back(ToRefine[i]);
			if(!directions.empty())
				newdirections.push_back(directions[i]);
		}
	}
	directions=newdirections;
	return res;
}

vector<int> RefineStrategy::RemoveDuplicatedLately(vector<int> const& ToRefine,
	int Npoints,vector<Vector2D> &directions,vector<int> const& Removed,
	Tessellation const& tess)
{
	if(!Removed.empty())
	{
		// Update the refined_old list
		size_t toAdd;
		for(size_t i=0;i<refined_old.size();++i)
		{
			toAdd=lower_bound(Removed.begin(),Removed.end(),refined_old[i])-
				Removed.begin();
			refined_old[i]-=int(toAdd);
		}
	}
	vector<int> res;
	vector<Vector2D> newdirections;
	for(size_t i=0;i<ToRefine.size();++i)
	{
		if(!binary_search(refined_old.begin(),refined_old.end(),ToRefine[i]))
		{
			// Make sure not to close to any neighbor
			const double R=tess.GetWidth(ToRefine[i]);
			vector<int> neigh=tess.GetNeighbors(ToRefine[i]);
			RemoveVal(neigh,-1);
			const int nn=(int)neigh.size();
			bool good=true;
			for(int j=0;j<nn;++j)
				if(tess.GetMeshPoint(neigh[j]).distance(
					tess.GetMeshPoint(ToRefine[i]))<0.2*R)
					good=false;
			if(good)
			{
				res.push_back(ToRefine[i]);
				if(!directions.empty())
					newdirections.push_back(directions[i]);
			}
		}
	}
	directions=newdirections;
	vector<int> temp=res;
	sort(temp.begin(),temp.end());
	refined_old=temp;
	// Add the new points that will be added
	int N=int(res.size());
	for(int i=0;i<N;++i)
		refined_old.push_back(Npoints+i);
	return res;
}

RefineStrategy::~RefineStrategy(void) {}

Vector2D FindBestSplit(Tessellation const* tess,int PointToRefine,
	vector<Edge> const& edges,double R,Vector2D &normal)
{
	Vector2D slope;
	Vector2D point(tess->GetMeshPoint(PointToRefine));
	int nedges=(int) edges.size();
	if(point.distance(tess->GetCellCM(PointToRefine))<0.1*R)
	{
		// Do Nearest edge
		double dis=DistanceToEdge(point,edges[0]);
		int min_edge=0;
		for(int j=1;j<nedges;++j)
		{
			double temp=DistanceToEdge(point,edges[j]);
			if(temp<dis)
			{
				dis=temp;
				min_edge=j;
			}
		}
		slope=Parallel(edges[min_edge]);
		slope=slope/abs(slope);
		Vector2D v=point-edges[min_edge].vertices.first;
		normal=v-ScalarProd(slope,v)*slope;
		normal=normal/abs(normal);
	}
	else
	{
		normal=tess->GetCellCM(PointToRefine)-point;
		normal=normal/abs(normal);
		slope.Set(normal.y,-normal.x);
	}
	return slope;
}
