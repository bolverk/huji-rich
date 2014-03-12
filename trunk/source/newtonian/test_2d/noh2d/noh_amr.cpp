#include "noh_amr.hpp"

NohRefine::NohRefine(double Vmax):Vmax_(Vmax)
{}

NohRefine::~NohRefine(void)
{}

vector<int> NohRefine::CellsToRefine(Tessellation const& tess,
	vector<Primitive> const& cells,vector<vector<double> > const& /*tracers*/,
	double /*time*/,vector<Vector2D> &directions,vector<int> const& Removed)
{
	vector<int> res;
	directions.clear();
	for(int i=0;i<(int)cells.size();++i)
	{
		// Are we too big? If so add to the refinement list
		if(tess.GetVolume(i)>Vmax_)
		{
			res.push_back(i);
			// Give a prefered splitting direction, this is not really needed in this case since the default scheme works very well
			Vector2D point=tess.GetMeshPoint(i);
			directions.push_back(Vector2D(point.x,point.y)/abs(point));
		}
	}
	// Make sure we are not splitting a cell turn after turn
	return RemoveDuplicatedLately(res,tess.GetPointNo(),directions,
		Removed);
}

NohRemove::NohRemove(double Vmin):Vmin_(Vmin)
{}

NohRemove::~NohRemove()
{}

vector<int> NohRemove::CellsToRemove(Tessellation const& tess,
	vector<Primitive> const& cells,vector<vector<double> > const& /*tracers*/,
	double /*time*/)const
{
	vector<int> ToRemoveTemp;
	vector<double> merit;
	for(int i=0;i<(int)cells.size();++i)
	{
		if(tess.GetVolume(i)<Vmin_)
		{
			ToRemoveTemp.push_back(i);
			merit.push_back(tess.GetVolume(i));
		}
	}
	// Make sure there are no neighbors
	vector<int> ToRemove=RemoveNeighbors(merit,ToRemoveTemp,tess);
	// Be 100% sure we removed the neighbors correctly
	CheckOutput(tess,ToRemove);
	return ToRemove;
}
