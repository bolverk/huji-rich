#include "RemovalStrategy.hpp"

vector<int> RemovalStrategy::RemoveNeighbors
(vector<double> const& merits,vector<int> const& 
 candidates,Tessellation const& tess) const
{
	vector<int> result;
	if(merits.size()!=candidates.size())
		throw UniversalError("Merits and Candidates don't have same size in RemovalStrategy");
	// Make sure there are no neighbors
	vector<int> bad_neigh;
	  int n=(int)merits.size();
	for(int i=0;i<n;++i)
	{
		vector<int> neigh=tess.GetNeighbors(candidates[i]);
		bool good=true;
		if(find(bad_neigh.begin(),bad_neigh.end(),candidates[i])!=
			bad_neigh.end())
			good=false;
		else
		{
			for(int j=0;j<(int)neigh.size();++j)
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
			result.push_back(candidates[i]);
	}
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
/*			if(temp.GetNeighbor(0)>N||temp.GetNeighbor(1)>N)
				throw UniversalError("Bad removal, neighbor is periodic");*/
			if(temp.GetNeighbor(0)==ToRemove[i])
			{
				if(binary_search(ToRemove.begin(),ToRemove.end(),temp.
					GetNeighbor(1)))
					throw UniversalError("Bad removal, neighboring cells");
			}
			else
				if(binary_search(ToRemove.begin(),ToRemove.end(),temp.
					GetNeighbor(0)))
					throw UniversalError("Bad removal, neighboring cells");
		}
	}
	return;
}

RemovalStrategy::~RemovalStrategy(void) {}
