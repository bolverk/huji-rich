#include "ConvexHull.hpp"

namespace
{
	bool VectorSort(Vector2D const& v1,Vector2D const& v2)
	{
		return ((v1.y<v2.y)||((v1.y==v2.y)&&(v1.x<v2.x)));
	}

	bool AngleSort(Vector2D const& v1,Vector2D const& v2)
	{
		const double tol=1e-8;
		const double a1=atan2(v1.y,v1.x);
		const double a2=atan2(v2.y,v2.x);
		if(a1<a2+tol)
			return true;
		if(a2<a1+tol)
			return false;
		if(abs(v1)<abs(v2))
			return true;
		else
			return false;
	}

	Vector2D GetReflectedPoint(Edge const& edge,Vector2D const& point)
	{
		Vector2D par(Parallel(edge));
		par=par/abs(par);
		Vector2D edge0=edge.GetVertex(0);
		Vector2D temp=point-edge0;
		return 2*par*ScalarProd(par,temp)-temp+edge0;
	}
}


void ConvexHull(vector<Vector2D> &result,Tessellation const* tess,int index)
{
	vector<int> edge_index=tess->GetCellEdges(index);
	const double eps=1e-7;
	vector<Vector2D> points;
	double R=tess->GetWidth(index);
	points.push_back(tess->GetEdge(edge_index[0]).GetVertex(0));
	points.push_back(tess->GetEdge(edge_index[0]).GetVertex(1));
	// Remove identical points
	for(size_t i=1;i<edge_index.size();++i)
	{
		size_t n=points.size();
		bool samepoint=false;
		for(size_t j=0;j<n;++j)
		{
			if(tess->GetEdge(edge_index[i]).GetVertex(0).distance(points[j])<eps*R)
				samepoint=true;
		}
		if(!samepoint)
			points.push_back(tess->GetEdge(edge_index[i]).GetVertex(0));
		samepoint=false;
		for(size_t j=0;j<n;++j)
		{
			if(tess->GetEdge(edge_index[i]).GetVertex(1).distance(points[j])<eps*R)
				samepoint=true;
		}
		if(!samepoint)
			points.push_back(tess->GetEdge(edge_index[i]).GetVertex(1));
	}
	// Find the bottom point
	sort(points.begin(),points.end(),VectorSort);
	// Start building the convexhull
	int n=(int)points.size();
	vector<int> indeces(n-1);
	vector<double> angles(n-1);
	for(int i=0;i<n-1;++i)
		angles[i]=atan2(points[i+1].y-points[0].y,points[i+1].x-points[0].x);
	sort_index(angles,indeces);
	result.resize(points.size());
	result[0]=points[0];
	// Check for colinear points
	const double tol=1e-8;
	vector<Vector2D> pfirst,plast;
	for(int i=1;i<n-1;++i)
	{
		if(abs(angles[indeces[i]]-angles[indeces[0]])<tol)
			pfirst.push_back(points[indeces[i]+1]);
		if(abs(angles[indeces[n-i-2]]-angles[indeces[n-2]])<tol)
			plast.push_back(points[indeces[n-i-2]+1]);
	}
	result[0]=points[0];
	if(!pfirst.empty())
	{
		pfirst.insert(pfirst.begin(),points[indeces[0]+1]);
		int N=(int)pfirst.size();
		vector<double> dist(N);
		for(int i=0;i<N;++i)
			dist[i]=abs(pfirst[i]-points[0]);
		vector<int> indeces2(N);
		sort_index(dist,indeces2);
		ReArrangeVector(pfirst,indeces2);
		for(int i=0;i<N;++i)
			result[i+1]=pfirst[i];
	}
	if(!plast.empty())
	{
		plast.insert(plast.begin(),points[indeces[n-2]+1]);
		int N=(int)plast.size();
		vector<double> dist;
		for(int i=0;i<N;++i)
			dist.push_back(abs(plast[i]-points[0]));
		vector<int> indeces2(N);
		sort_index(dist,indeces2);
		ReArrangeVector(plast,indeces2);
		for(int i=0;i<N;++i)
			result[n-1-i]=plast[i];
	}
	int loc1=(int)pfirst.size();
	int loc2=(int)plast.size();
	for(int i=loc1+1;i<n-loc2;++i)
		result[i]=points[indeces[i-1]+1];
}

void ConvexEdges(vector<int> &result,Tessellation const* tess,int index)
{
	vector<int> const& edges=tess->GetCellEdges(index);
	const Vector2D mypoint=tess->GetMeshPoint(index);
	int nedges=(int)edges.size();
	result.resize(nedges);
	vector<double> angles(nedges);
	for(int i=0;i<nedges;++i)
	{
		Edge const& edge=tess->GetEdge(edges[i]);
		const int other=(edge.GetNeighbor(0)==index)? edge.GetNeighbor(1):edge.GetNeighbor(0);
		Vector2D otherpoint=(other==-1) ? GetReflectedPoint(edge,mypoint) : tess->GetMeshPoint(other);
		angles[i]=atan2(otherpoint.y-mypoint.y,otherpoint.x-mypoint.x);
	}
	vector<int> temp;
	sort_index(angles,temp);
	for(int i=0;i<nedges;++i)
		result[i]=edges[temp[i]];
}