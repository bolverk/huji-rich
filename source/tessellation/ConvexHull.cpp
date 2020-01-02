#include <boost/foreach.hpp>
#include <cfloat>
#include "ConvexHull.hpp"
#include "geotests.hpp"

using namespace std;

namespace
{
	Vector2D GetReflectedPoint(Edge const& edge,Vector2D const& point)
	{
	  const Vector2D par(normalize(Parallel(edge)));
	  const Vector2D edge0=edge.vertices.first;
	  const Vector2D temp=point-edge0;
	  return 2*par*ScalarProd(par,temp)-temp+edge0;
	}

  bool check_same_point(const vector<Vector2D>& vertices,
			const Vector2D& p,
			double tol)
  {
    BOOST_FOREACH(const Vector2D& v, vertices)
      {
	if(dist_sqr(v-p)<tol)
	  return true;
      }
    return false;
  }
}

void ConvexHull(vector<Vector2D> &result,Tessellation const& tess,int index)
{
	vector<int> edge_index=tess.GetCellEdges(index);
	const double eps=1e-30;
	vector<Vector2D> points;
	points.reserve(10);
	double R=tess.GetWidth(index);
	points.push_back(tess.GetEdge(edge_index[0]).vertices.first);
	points.push_back(tess.GetEdge(edge_index[0]).vertices.second);
	// Remove bad edges
	for(size_t i=1;i<edge_index.size();++i)
	{
	  const Edge& edge = tess.GetEdge(edge_index[i]);
	  if (edge.GetLength() < 1e-8 * R)
		  continue;
		if(!check_same_point(points,edge.vertices.first,eps*pow(R,2)))
		  points.push_back(edge.vertices.first);
		if(!check_same_point(points,edge.vertices.second,eps*pow(R,2)))
		  points.push_back(edge.vertices.second);
	}

	const Vector2D cm = tess.GetCellCM(index);
	
	// Start building the convexhull
	size_t n=points.size();
	vector<double> angles(n);
	for(size_t i=0;i<n;++i)
	  angles.at(i)=atan2(points.at(i).y-cm.y,points.at(i).x-cm.x);
	const vector<size_t> indeces = sort_index(angles);
	result = VectorValues(points,indeces);
}

void ConvexHull(std::vector<std::pair<Vector2D, Vector2D> >& result, Tessellation const& tess, int index)
{
	vector<int> edge_index = tess.GetCellEdges(index);
	size_t Nedges = edge_index.size();
	double R = 0.0;
	double eps = 1e-8;
	for (size_t i = 0; i < Nedges; ++i)
		R = std::max(R, tess.GetEdge(edge_index[i]).GetLength());
	// keep only good edges
	std::vector<int> good_edges;
	for (size_t i = 0; i < Nedges; ++i)
		if (tess.GetEdge(edge_index[i]).GetLength() > R* eps)
			good_edges.push_back(edge_index[i]);
	Nedges = good_edges.size();
	std::vector<double> angles(Nedges);
	const Vector2D cm = tess.GetCellCM(index);
	// sort by angle
	for (size_t i = 0; i < Nedges; ++i)
	{
		Vector2D midpoint = 0.5 * (tess.GetEdge(good_edges[i]).vertices.first + tess.GetEdge(good_edges[i]).vertices.second) - cm;
		angles[i] = atan2(midpoint.y, midpoint.x);
	}
	const vector<size_t> indeces = sort_index(angles);
	result.resize(Nedges);
	for (size_t i = 0; i < Nedges; ++i)
	{
		Edge const& edge = tess.GetEdge(good_edges[indeces[i]]);
		// make sure in right order
		if(orient2d(TripleConstRef<Vector2D>(edge.vertices.first, edge.vertices.second, cm)) > 0.0)
			result[i] = std::pair<Vector2D, Vector2D>(edge.vertices.first, edge.vertices.second);
		else
			result[i] = std::pair<Vector2D, Vector2D>(edge.vertices.second, edge.vertices.first);
	}
}

void ConvexEdges(vector<int> &result,Tessellation const& tess,int index)
{
	vector<int> const& edges=tess.GetCellEdges(index);
	const Vector2D mypoint=tess.GetMeshPoint(index);
	int nedges=static_cast<int>(edges.size());
	result.resize(static_cast<size_t>(nedges));
	vector<double> angles(static_cast<size_t>(nedges));
	for(int i=0;i<nedges;++i)
	{
		Edge const& edge=tess.GetEdge(edges[static_cast<size_t>(i)]);
		const int other=(edge.neighbors.first==index)? edge.neighbors.second : edge.neighbors.first;
		Vector2D otherpoint=(other==-1) ? GetReflectedPoint(edge,mypoint) : tess.GetMeshPoint(other);
		angles[static_cast<size_t>(i)]=atan2(otherpoint.y-mypoint.y,otherpoint.x-mypoint.x);
	}
	vector<int> temp;
	sort_index(angles,temp);
	for(size_t i=0;i<static_cast<size_t>(nedges);++i)
	  result[i]=edges[static_cast<size_t>(temp[i])];
}
