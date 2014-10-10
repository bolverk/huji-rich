#include "RemovalStrategy.hpp"

vector<int> RemovalStrategy::RemoveNearBoundary(vector<int> const& ToRemove,Tessellation
	const& tess)const
{
	int nrefine=(int)ToRemove.size();
	int npoints=tess.GetPointNo();
	vector<int> res;
	for(int i=0;i<nrefine;++i)
	{
		vector<int> const& edges=tess.GetCellEdges(ToRemove[i]);
		int nedge=(int) edges.size();
		bool good=true;
		for(int j=0;j<nedge;++j)
		{
			Edge const& edge=tess.GetEdge(edges[j]);
			if(edge.neighbors.first>npoints||edge.neighbors.second>npoints)
			{
				good=false;
				break;
			}
		}
		if(good)
			res.push_back(ToRemove[i]);
	}
	return res;
}

vector<int> RemovalStrategy::RemoveNeighbors
(vector<double> const& merits,vector<int> const&
 candidates,Tessellation const& tess) const
{
	vector<int> result;
	vector<double> merits2;
	if(merits.size()!=candidates.size())
		throw UniversalError("Merits and Candidates don't have same size in RemovalStrategy");
	// Make sure there are no neighbors
	vector<int> bad_neigh;
	int n=(int)merits.size();
	//	int npoints=tess.GetPointNo();
	for(int i=0;i<n;++i)
	{
		bool good=true;
		vector<int> edges=tess.GetCellEdges(candidates[i]);
		vector<int> neigh=tess.GetNeighbors(candidates[i]);
		int nneigh=(int) neigh.size();
		/*for(int j=0;j<nneigh;++j)
		{
			Edge const& edge=tess.GetEdge(edges[j]);
			if(edge.neighbors.first>npoints||edge.neighbors.second>npoints)
			{
				good=false;
				break;
			}
		}*/
		if(!good)
			continue;
		if(find(bad_neigh.begin(),bad_neigh.end(),candidates[i])!=
			bad_neigh.end())
			good=false;
		else
		{
			for(int j=0;j<nneigh;++j)
			{
				if(binary_search(candidates.begin(),candidates.end(),neigh[j]))
				{
					if(merits[i]<merits[lower_bound(candidates.begin(),
						candidates.end(),neigh[j])-candidates.begin()])
					{
						good=false;
						break;
					}
					if(merits[i]==merits[lower_bound(candidates.begin(),
						candidates.end(),neigh[j])-candidates.begin()])
					{
						if(find(bad_neigh.begin(),bad_neigh.end(),neigh[j])==
							bad_neigh.end())
							bad_neigh.push_back(neigh[j]);
					}
				}
			}
		}
		if(good)
		{
			result.push_back(candidates[i]);
			merits2.push_back(merits[i]);
		}
	}
#ifdef RICH_MPI
	result=RemoveMPINeighbors(result,merits2,tess);
#endif
	return result;
}

void RemovalStrategy::CheckOutput(Tessellation const& tess,vector<int>
	& ToRemove)const
{
	sort(ToRemove.begin(),ToRemove.end());
	int n=int(ToRemove.size());
	for(int i=0;i<n;++i)
	{
		vector<int> edges=tess.GetCellEdges(ToRemove[i]);
		//check we are not near periodic boundary
		for(int j=0;j<(int)edges.size();++j)
		{
			Edge temp=tess.GetEdge(edges[j]);
/*			if(temp.neighbors.first>N||temp.neighbors.second>N)
				throw UniversalError("Bad removal, neighbor is periodic");*/
			if(temp.neighbors.first==ToRemove[i])
			{
				if(binary_search(ToRemove.begin(),ToRemove.end(),temp.
					neighbors.second))
					throw UniversalError("Bad removal, neighboring cells");
			}
			else
				if(binary_search(ToRemove.begin(),ToRemove.end(),temp.
					neighbors.first))
					throw UniversalError("Bad removal, neighboring cells");
		}
	}
	return;
}

RemovalStrategy::~RemovalStrategy(void) {}
