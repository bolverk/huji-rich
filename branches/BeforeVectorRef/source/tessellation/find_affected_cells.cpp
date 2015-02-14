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

  /*
	vector<int> join(const vector<int>& v,
		int n)
	{
		vector<int> res = v;
		res.push_back(n);
		return res;
	}
  */

	double bracket(double low, double num, double high)
	{
		return std::max(std::min(num,high),low);
	}

	bool edge_circle_intersect(const Edge& edge,
		const Circle& circle)
	{
		const double s = ScalarProd(circle.getCenter() - edge.vertices.first,
			edge.vertices.second - edge.vertices.first)/
			ScalarProd(edge.vertices.second - edge.vertices.first,
			edge.vertices.second - edge.vertices.first);
		const double sb = bracket(0,s,1);
		return circle((1-sb)*edge.vertices.first+sb*edge.vertices.second);
	}

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
			BOOST_FOREACH(int nbr, tess.GetNeighbors(index))
			{			
				if(nbr!=-1)
					find_affected_cells2(tess,nbr,circle,res,visited);
				else
				{
					if(!contains(visited,-1))
						res.push_back(-1);
				}
			}
			return;
		}
	}
}

void find_affected_cells(const Tessellation& tess,int index,const Circle& circle,
	vector<int> &res)
{
	vector<int> visited;
	find_affected_cells2(tess,index,circle,res,visited);
}
