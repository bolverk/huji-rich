#include "EdgeLengthCorrect.hpp"

void CorrectEdgeLength(Tessellation const& tessold,Tessellation const& tessnew,
	vector<double> &lengths)
{
	int n=tessold.GetTotalSidesNumber();
	int npoints=tessold.GetPointNo();
	lengths.resize(static_cast<size_t>(n));
	for(int i=0;i<npoints;++i)
	{
		vector<int> edgesold=tessold.GetCellEdges(i);
		vector<int> edgesnew=tessnew.GetCellEdges(i);
		int nedges=static_cast<int>(edgesold.size());
		int nedgesnew=static_cast<int>(edgesnew.size());
		for(int j=0;j<nedges;++j)
		{
		  Edge const& edge=tessold.GetEdge(edgesold[static_cast<size_t>(j)]);
			int n0=edge.neighbors.first;
			int n1=edge.neighbors.second;
			bool found=false;
			// Don't do double work
			if((n0<i&&n0>=0)||(n1<i&&n1>=0))
				continue;
			for(int k=0;k<nedgesnew;++k)
			{
			  Edge const& edgenew=tessnew.GetEdge(edgesnew[static_cast<size_t>((k+j)%nedgesnew)]);
				// are the two edges the same?
				if(((tessold.GetOriginalIndex(n0)==tessnew.GetOriginalIndex(
					edgenew.neighbors.first))&&
					(tessold.GetOriginalIndex(n1)==tessnew.GetOriginalIndex(
					edgenew.neighbors.second)))||((tessold.GetOriginalIndex(n0)==
					tessnew.GetOriginalIndex(edgenew.neighbors.second))&&
					(tessold.GetOriginalIndex(n1)==tessnew.GetOriginalIndex(
					edgenew.neighbors.first))))
				{
				  lengths[static_cast<size_t>(edgesold[static_cast<size_t>(j)])]=0.5*(edge.GetLength()+edgenew.GetLength());
					found=true;
					break;
				}
			}
			if(!found)
			  lengths[static_cast<size_t>(edgesold[static_cast<size_t>(j)])]=0.5*edge.GetLength();
		}
	}
}

void CorrectEdgeLength(Tessellation const& tessold,Tessellation const& tessmid,
	Tessellation const& tessnew,vector<double> &lengths)
{
	int n=tessmid.GetTotalSidesNumber();
	int npoints=tessmid.GetPointNo();
	lengths.resize(static_cast<size_t>(n));
	for(int i=0;i<npoints;++i)
	{
		vector<int> edgesold=tessold.GetCellEdges(i);
		vector<int> edgesnew=tessnew.GetCellEdges(i);
		vector<int> edgesmid=tessmid.GetCellEdges(i);
		int nedges=static_cast<int>(edgesold.size());
		int nedgesnew=static_cast<int>(edgesnew.size());
		int nedgesmid=static_cast<int>(edgesmid.size());
		for(int j=0;j<nedgesmid;++j)
		{
		  Edge const& edge=tessmid.GetEdge(edgesmid[static_cast<size_t>(j)]);
			int n0=edge.neighbors.first;
			int n1=edge.neighbors.second;
			// Don't do double work
			if((n0<i&&n0>=0)||(n1<i&&n1>=0))
				continue;
			lengths[static_cast<size_t>(edgesmid[static_cast<size_t>(j)])]=0.5*edge.GetLength();
			for(int k=0;k<nedgesnew;++k)
			{
			  Edge const& edgenew=tessnew.GetEdge(edgesnew[static_cast<size_t>((k+j)%nedgesnew)]);
				// are the two edges the same?
				if(((tessmid.GetOriginalIndex(n0)==tessnew.GetOriginalIndex(
					edgenew.neighbors.first))&&
					(tessmid.GetOriginalIndex(n1)==tessnew.GetOriginalIndex(
					edgenew.neighbors.second)))||((tessmid.GetOriginalIndex(n0)==
					tessnew.GetOriginalIndex(edgenew.neighbors.second))&&
					(tessmid.GetOriginalIndex(n1)==tessnew.GetOriginalIndex(
					edgenew.neighbors.first))))
				{
				  lengths[static_cast<size_t>(edgesmid[static_cast<size_t>(j)])]+=0.25*edgenew.GetLength();
					break;
				}
			}
			for(int k=0;k<nedges;++k)
			{
			  Edge const& edgeold=tessold.GetEdge(edgesold[static_cast<size_t>((k+j)%nedges)]);
				// are the two edges the same?
				if(((tessmid.GetOriginalIndex(n0)==tessold.GetOriginalIndex(
					edgeold.neighbors.first))&&
					(tessmid.GetOriginalIndex(n1)==tessold.GetOriginalIndex(
					edgeold.neighbors.second)))||((tessmid.GetOriginalIndex(n0)==
					tessold.GetOriginalIndex(edgeold.neighbors.second))&&
					(tessmid.GetOriginalIndex(n1)==tessold.GetOriginalIndex(
					edgeold.neighbors.first))))
				{
				  lengths[static_cast<size_t>(edgesmid[static_cast<size_t>(j)])]+=0.25*edgeold.GetLength();
					break;
				}
			}
		}
	}
}
