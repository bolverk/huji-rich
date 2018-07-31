#include "find_affected_cells.hpp"
#include "geometry.hpp"
#include "boost/foreach.hpp"
#include "../misc/utils.hpp"

namespace 
{
	bool periodic_duplicate_before(std::vector<Vector2D> const& toadd,Vector2D const& tocheck)
	{
		size_t N = toadd.size();
		if (N == 0)
			return false;
		double dx_1 = 0;
		for (size_t i = 0; i < N; ++i)
			dx_1 = std::max(dx_1, 1.0 / fastabs(toadd[i]));
		double dx = 1.0 / dx_1;
		for (size_t i = 0; i < N; ++i)
		{
			double dV = fastabs(toadd[i] - tocheck);
			if (dV < dx*0.1)
			{
				return true;
			}
		}
		return false;
	}

	bool contains(const vector<int>& v,
		const int n)
	{
		for(size_t i=0;i<v.size();++i){
			if(n==v[i])
				return true;
		}
		return false;
	}

	double bracket(double low, double num, double high)
	{
		return std::max(std::min(num,high),low);
	}
}

bool edge_circle_intersect
(const Edge& edge,
 const Circle& circle)
{
  const double s = 
    ScalarProd(circle.getCenter() - edge.vertices.first,
	       edge.vertices.second - edge.vertices.first)/
    ScalarProd(edge.vertices.second - edge.vertices.first,
	       edge.vertices.second - edge.vertices.first);
  const double sb = bracket(0,s,1);
  return circle((1-sb)*edge.vertices.first+sb*edge.vertices.second);
}

namespace {
	bool cell_circle_intersect(const Tessellation& tess,
		int index,
		const Circle& circle)
	{
		BOOST_FOREACH(int ei, tess.GetCellEdges(index))
		{
			if(edge_circle_intersect(tess.GetEdge(ei),circle))
				return true;
		}
		return false;
	}

	void find_affected_cells2(const Tessellation& tess,int index, Circle circle,
		vector<int> &res,vector<int>& visited,std::vector<Vector2D> &added,bool periodic,Vector2D const& toaddsingle)
	{
		if(contains(visited,index))
			return;
		visited.push_back(index);
		if(!cell_circle_intersect(tess,tess.GetOriginalIndex(index),circle))
			return;
		else
		{
			res.push_back(index);
			if (periodic)
				added.push_back(toaddsingle);
			vector<int> vtemp = tess.GetNeighbors(index);
			size_t NN = vtemp.size();
			for(size_t j = 0; j < NN;++j)
			{			
			  if(vtemp[j]<tess.GetPointNo())
					find_affected_cells2(tess, vtemp[j],circle,res,visited,added,periodic,toaddsingle);
			  else
				  if (periodic)
				  {
					  const Vector2D toadd = tess.GetMeshPoint(tess.GetOriginalIndex(vtemp[j])) -
						  tess.GetMeshPoint(vtemp[j]);
					  const Vector2D toaddsinglenew = toaddsingle + toadd;
					  // Make sure not heading back in duplicate space
					  if (!periodic_duplicate_before(added, toadd))
					  {
						  Vector2D oldcenter = circle.getCenter();
						  circle.setCenter(toadd + oldcenter);						  
						  find_affected_cells2(tess, tess.GetOriginalIndex(vtemp[j]), circle, res, visited, added, periodic,
							  toaddsinglenew);
						  circle.setCenter(oldcenter);
					  }			  
				  }
			}
			return;
		}
	}
}

vector<int> find_affected_cells(const Tessellation& tess, int index, Circle& circle, vector<int> &vtemp,bool periodic,
	std::vector<Vector2D> &periodic_add)
{
	periodic_add.clear();
  vector<int> res;
  vtemp = tess.GetNeighbors(index);
  size_t N = vtemp.size();
  for (size_t i = 0; i < N;++i)
  {
	  
		  if (vtemp[i] < tess.GetPointNo())
		  {
			  if (cell_circle_intersect(tess, vtemp[i], circle))
			  {
				  res.push_back(vtemp[i]);
				  if (periodic)
					  periodic_add.push_back(Vector2D());
			  }
		  }
		  else
		  {
			  if (periodic)
			  {
				  int real_neigh = tess.GetOriginalIndex(vtemp[i]);
				  Vector2D toadd = tess.GetMeshPoint(vtemp[i]) - tess.GetMeshPoint(real_neigh);
				  circle.setCenter(circle.getCenter() + toadd);
				  if (cell_circle_intersect(tess, real_neigh, circle))
				  {
					  res.push_back(real_neigh);
					  periodic_add.push_back(toadd);
				  }
				  circle.setCenter(circle.getCenter() - toadd);
			  }
		  }
  }
  return res;
}

void find_affected_cells_recursive(const Tessellation& tess,int index,const Circle& circle,
	vector<int> &res,std::vector<Vector2D> &added,bool periodic)
{
	vector<int> visited;
	added.clear();
	Vector2D toaddsingle;
	find_affected_cells2(tess,index,circle,res,visited,added,periodic,toaddsingle);
}
