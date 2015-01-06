#include <cfloat>
#include "ConvexHull.hpp"

using namespace std;

namespace
{
	bool VectorSort(Vector2D const& v1,Vector2D const& v2)
	{
		const double relval=max(max(fabs(v1.y),fabs(v2.y)),DBL_EPSILON*2);
		const double diffmax=4*max(fabs(v1.y-v2.y),fabs(v2.y-v1.y))/relval;
		if(diffmax<DBL_EPSILON)
			return v1.x<v2.x;
		else
			return v1.y<v2.y;
	}

	/*
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
	*/

	Vector2D GetReflectedPoint(Edge const& edge,Vector2D const& point)
	{
		Vector2D par(Parallel(edge));
		par=par/abs(par);
		Vector2D edge0=edge.vertices.first;
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
	points.push_back(tess->GetEdge(edge_index[0]).vertices.first);
	points.push_back(tess->GetEdge(edge_index[0]).vertices.second);
	// Remove identical points
	for(size_t i=1;i<edge_index.size();++i)
	{
		size_t n=points.size();
		bool samepoint=false;
		for(size_t j=0;j<n;++j)
		{
			if(tess->GetEdge(edge_index[i]).vertices.first.distance(points[j])<eps*R)
				samepoint=true;
		}
		if(!samepoint)
			points.push_back(tess->GetEdge(edge_index[i]).vertices.first);
		samepoint=false;
		for(size_t j=0;j<n;++j)
		{
			if(tess->GetEdge(edge_index[i]).vertices.second.distance(points[j])<eps*R)
				samepoint=true;
		}
		if(!samepoint)
			points.push_back(tess->GetEdge(edge_index[i]).vertices.second);
	}
	// Find the bottom point
	sort(points.begin(),points.end(),VectorSort);
	// Start building the convexhull
	size_t n=points.size();
	//	vector<int> indeces(n-1);
	vector<double> angles(n-1);
	for(size_t i=0;i<n-1;++i)
		angles[i]=atan2(points[i+1].y-points[0].y,points[i+1].x-points[0].x);
	//	sort_index(angles,indeces);
	const vector<size_t> indeces = sort_index(angles);
	result.resize(points.size());
	result[0]=points[0];
	// Check for colinear points
	const double tol=1e-8;
	vector<Vector2D> pfirst,plast;
	for(size_t i=1;i<n-1;++i)
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
		//		int N=(int)pfirst.size();
		vector<double> dist(pfirst.size());
		for(size_t i=0;i<pfirst.size();++i)
			dist[i]=abs(pfirst[i]-points[0]);
		vector<int> indeces2(pfirst.size());
		sort_index(dist,indeces2);
		ReArrangeVector(pfirst,indeces2);
		for(size_t i=0;i<pfirst.size();++i)
			result[i+1]=pfirst[i];
	}
	if(!plast.empty())
	{
		plast.insert(plast.begin(),points[(size_t)indeces[(size_t)n-2]+1]);
		int N=(int)plast.size();
		vector<double> dist;
		for(int i=0;i<N;++i)
			dist.push_back(abs(plast[(size_t)i]-points[0]));
		vector<int> indeces2((size_t)N);
		sort_index(dist,indeces2);
		ReArrangeVector(plast,indeces2);
		for(int i=0;i<N;++i)
			result[(size_t)(n-1-i)]=plast[(size_t)i];
	}
	int loc1=(int)pfirst.size();
	int loc2=(int)plast.size();
	for(size_t i=loc1+1;i<n-loc2;++i)
	  result[i]=points[static_cast<size_t>(indeces[i-1])+1];
}

void ConvexEdges(vector<int> &result,Tessellation const* tess,int index)
{
	vector<int> const& edges=tess->GetCellEdges(index);
	const Vector2D mypoint=tess->GetMeshPoint(index);
	int nedges=(int)edges.size();
	result.resize((size_t)nedges);
	vector<double> angles((size_t)nedges);
	for(int i=0;i<nedges;++i)
	{
		Edge const& edge=tess->GetEdge(edges[(size_t)i]);
		const int other=(edge.neighbors.first==index)? edge.neighbors.second : edge.neighbors.first;
		Vector2D otherpoint=(other==-1) ? GetReflectedPoint(edge,mypoint) : tess->GetMeshPoint(other);
		angles[(size_t)i]=atan2(otherpoint.y-mypoint.y,otherpoint.x-mypoint.x);
	}
	vector<int> temp;
	sort_index(angles,temp);
	for(size_t i=0;i<(size_t)nedges;++i)
		result[i]=edges[(size_t)temp[i]];
}
