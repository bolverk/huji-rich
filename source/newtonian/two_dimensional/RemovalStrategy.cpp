#include "RemovalStrategy.hpp"
#include <cmath>

vector<int> RemovalStrategy::RemoveNearBoundary(vector<int> const& ToRemove,Tessellation
	const& tess)const
{
	int nrefine=static_cast<int>(ToRemove.size());
	int npoints=tess.GetPointNo();
	vector<int> res;
	for(int i=0;i<nrefine;++i)
	{
	  vector<int> const& edges=tess.GetCellEdges(ToRemove[static_cast<size_t>(i)]);
	  int nedge=static_cast<int>(edges.size());
		bool good=true;
		for(int j=0;j<nedge;++j)
		{
		  Edge const& edge=tess.GetEdge(edges[static_cast<size_t>(j)]);
			if(edge.neighbors.first>npoints||edge.neighbors.second>npoints)
			{
				good=false;
				break;
			}
		}
		if(good)
		  res.push_back(ToRemove[static_cast<size_t>(i)]);
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
	int n=static_cast<int>(merits.size());
	//	int npoints=tess.GetPointNo();
	for(int i=0;i<n;++i)
	{
		bool good=true;
		vector<int> neigh=tess.GetNeighbors(candidates[static_cast<size_t>(i)]);
		int nneigh=static_cast<int>(neigh.size());
		if(find(bad_neigh.begin(),bad_neigh.end(),candidates[static_cast<size_t>(i)])!=
			bad_neigh.end())
			good=false;
		else
		{
			for(int j=0;j<nneigh;++j)
			{
			  if(binary_search(candidates.begin(),candidates.end(),neigh[static_cast<size_t>(j)]))
				{
				  if(merits[static_cast<size_t>(i)]<merits[static_cast<size_t>(lower_bound(candidates.begin(),
										   candidates.end(),neigh[static_cast<size_t>(j)])-candidates.begin())])
					{
						good=false;
						break;
					}
				  if(fabs(merits[static_cast<size_t>(i)]-merits[static_cast<size_t>(lower_bound(candidates.begin(),
											candidates.end(),neigh[static_cast<size_t>(j)])-candidates.begin())])<1e-9)
					{
						if(find(bad_neigh.begin(),bad_neigh.end(),neigh[static_cast<size_t>(j)])==
							bad_neigh.end())
							bad_neigh.push_back(neigh[static_cast<size_t>(j)]);
					}
				}
			}
		}
		if(good)
		{
			result.push_back(candidates[static_cast<size_t>(i)]);
			merits2.push_back(merits[static_cast<size_t>(i)]);
		}
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
		vector<int> edges=tess.GetCellEdges(ToRemove[static_cast<size_t>(i)]);
		//check we are not near periodic boundary
		for(int j=0;j<static_cast<int>(edges.size());++j)
		{
			Edge temp=tess.GetEdge(edges[static_cast<size_t>(j)]);
/*			if(temp.neighbors.first>N||temp.neighbors.second>N)
				throw UniversalError("Bad removal, neighbor is periodic");*/
			if(temp.neighbors.first==ToRemove[static_cast<size_t>(i)])
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
