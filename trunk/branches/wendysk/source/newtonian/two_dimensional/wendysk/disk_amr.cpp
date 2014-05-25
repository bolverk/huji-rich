#include "disk_amr.hpp"

DiskRemove::DiskRemove(double inner_radius,
		       double outer_radius,
		       int total_specials):
  inner_radius_(inner_radius),
  outer_radius_(outer_radius),
  total_specials_(total_specials) {}

vector<int> DiskRemove::CellsToRemove
(Tessellation const* tess,
 vector<Primitive> const& /*cells*/,
 vector<vector<double> > const& /*tracers*/,
 double /*time*/) const
{
  vector<int> candidates;
  vector<double> merits;
  for(int i=total_specials_;i<tess->GetPointNo();++i){
    const Vector2D mp = tess->GetMeshPoint(i);
    if(abs(mp)<inner_radius_){
      candidates.push_back(i);
      merits.push_back(inner_radius_-abs(mp));
    }
    else if(abs(mp)>outer_radius_){
      candidates.push_back(i);
      merits.push_back(abs(mp)-outer_radius_);
    }
  }
  vector<int> result = RemoveNeighbors(merits,
				       candidates,
				       tess);
  CheckOutput(tess, result);
  return result;
}

DiskRefine::DiskRefine(int total_specials,
		       double d_min,
		       double min_radius,
			   double max_radius,
		       double max_volume):
  total_specials_(total_specials),
  d_min_(d_min),
  min_radius_(min_radius),
  max_radius_(max_radius),
  max_volume_(max_volume) {}

vector<int> DiskRefine::CellsToRefine
(Tessellation const* tess,
 vector<Primitive> const& cells,vector<vector<double> > const& /*tracers*/,
 double /*time*/,vector<Vector2D> & directions ,vector<int> const& Removed)
{
  vector<int> stage1;
  for(int i=total_specials_;i<tess->GetPointNo();++i){
    if(cells[i].Density>d_min_)
      stage1.push_back(i);
  }
  vector<int> stage2;
  for(int i=0;i<(int)stage1.size();++i){
    if((abs(tess->GetMeshPoint(stage1[i]))>min_radius_)&&
		abs(tess->GetMeshPoint(stage1[i]))<max_radius_)
      stage2.push_back(stage1[i]);
  }
   vector<int> res;
  for(int i=0;i<(int)stage2.size();++i){
    if(tess->GetVolume(stage2[i])>max_volume_
		*max(abs(tess->GetMeshPoint(stage2[i])),0.1)*3)
      res.push_back(stage2[i]);
  }

 
  return RemoveDuplicatedLately(res,
				tess->GetPointNo(),
				directions,
				Removed);
}
