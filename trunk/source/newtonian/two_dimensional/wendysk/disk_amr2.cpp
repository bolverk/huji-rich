#include "disk_amr2.hpp"

DiskRemove2::DiskRemove2(double inner_radius1,double inner_radius2,
	double outer_radius,int total_specials,Vector2D s1,Vector2D s2):
inner_radius1_(inner_radius1),inner_radius2_(inner_radius2),
	outer_radius_(outer_radius),
	total_specials_(total_specials),s1_(s1),s2_(s2) {}

vector<int> DiskRemove2::CellsToRemove(Tessellation const* tess,
	vector<Primitive> const& /*cells*/,vector<vector<double> > const& /*tracers*/,
	double /*time*/) const
{
	vector<int> candidates;
	vector<double> merits;
	for(int i=total_specials_;i<tess->GetPointNo();++i)
	{
		const Vector2D mp = tess->GetMeshPoint(i);
		if(abs(mp-s1_)<inner_radius1_)
		{
			candidates.push_back(i);
			merits.push_back(inner_radius1_-abs(mp));
			continue;
		}
		if(abs(mp)>outer_radius_)
		{
			candidates.push_back(i);
			merits.push_back(abs(mp)-outer_radius_);
			continue;
		}
		if(abs(mp-s2_)<inner_radius2_)
		{
			candidates.push_back(i);
			merits.push_back(inner_radius2_-abs(mp));
		}
	}
	vector<int> result = RemoveNeighbors(merits,
		candidates,
		tess);
	CheckOutput(tess, result);
	return result;
}

DiskRefine2::DiskRefine2(int total_specials,double d_min,double min_radius1,
	double min_radius2,double max_radius,double max_volume,Vector2D s1,
	Vector2D s2):
total_specials_(total_specials),d_min_(d_min),min_radius1_(min_radius1),
	min_radius2_(min_radius2),max_radius_(max_radius),
	max_volume_(max_volume),s1_(s1),s2_(s2) {}

vector<int> DiskRefine2::CellsToRefine(Tessellation const* tess,
	vector<Primitive> const& cells,vector<vector<double> > const& /*tracers*/,
	double /*time*/,vector<Vector2D> & directions ,vector<int> const& Removed)
{
	vector<int> stage1;
	for(int i=total_specials_;i<tess->GetPointNo();++i)
	{
		if(cells[i].Density>d_min_)
			stage1.push_back(i);
	}
	vector<int> stage2;
	for(int i=0;i<(int)stage1.size();++i)
	{
		if((abs(tess->GetMeshPoint(stage1[i])-s1_)>min_radius1_)&&
			abs(tess->GetMeshPoint(stage1[i]))<max_radius_)
		{
			stage2.push_back(stage1[i]);
			continue;
		}
		if((abs(tess->GetMeshPoint(stage1[i])-s2_)>min_radius2_)&&
			abs(tess->GetMeshPoint(stage1[i]))<max_radius_)
		{
			stage2.push_back(stage1[i]);
			continue;
		}
	}
	vector<int> res;
	for(int i=0;i<(int)stage2.size();++i)
	{
		if(tess->GetVolume(stage2[i])>max_volume_
			*max(min(abs(tess->GetMeshPoint(stage2[i])-s1_),
			abs(tess->GetMeshPoint(stage2[i])-s2_)),0.05)*3)
		{
			res.push_back(stage2[i]);
			continue;
		}
		if(tess->GetVolume(stage2[i])>(2*max_volume_))
			res.push_back(stage2[i]);
	}

	return RemoveDuplicatedLately(res,tess->GetPointNo(),directions,	Removed);
}
