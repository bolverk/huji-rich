#include "find_affected_cells.hpp"
#include "geometry.hpp"
#include "boost/foreach.hpp"
#include "../misc/utils.hpp"

namespace {
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

  /*
  int cell_circle_intersect_2
  (const Tessellation& tess,
   int index,
   const Circle& circle)
  {
    BOOST_FOREACH(int ei, tess.GetCellEdges(index)){
      if(edge_circle_intersect(tess.GetEdge(ei),circle))
	return ei;
    }
    return -1;
  }
  */

	void find_affected_cells2(const Tessellation& tess,int index,const Circle& circle,
		vector<int> &res,vector<int>& visited)
	{
		if(contains(visited,index))
			return;
		visited.push_back(index);
		if(!cell_circle_intersect(tess,index,circle))
			return;
		else
		{
			res.push_back(index);
			vector<int> vtemp = tess.GetNeighbors(index);
			size_t NN = vtemp.size();
			for(size_t j = 0; j < NN;++j)
			{			
			  if(vtemp[j]<tess.GetPointNo())
					find_affected_cells2(tess, vtemp[j],circle,res,visited);
			}
			return;
		}
	}
}

vector<int> find_affected_cells
(const Tessellation& tess,
 int index,
 const Circle& circle, vector<int> &vtemp)
{
  vector<int> res;
  vtemp = tess.GetNeighbors(index);
  size_t N = vtemp.size();
  for (size_t i = 0; i < N;++i)
  {
    if(vtemp[i]<tess.GetPointNo() && 
       cell_circle_intersect(tess,vtemp[i],circle))
      res.push_back(vtemp[i]);      
  }
  return res;
}

void find_affected_cells_recursive(const Tessellation& tess,int index,const Circle& circle,
	vector<int> &res)
{
	vector<int> visited;
	find_affected_cells2(tess,index,circle,res,visited);
}
