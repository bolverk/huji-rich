#include "SelfGravity2D.hpp"
#include "../../../misc/universal_error.hpp"

SelfGravity::SelfGravity(double OpenAngle, double hFactor):
  _counter(0),
  _length(0),
  _OpenAngle(OpenAngle),
  _tree(0),
  _treePoints(annAllocPts(2,2)),
  masses(vector<double>()),
  h(vector<double>()),
  _Factor(hFactor) {}

SelfGravity::~SelfGravity()
{
}

Vector2D SelfGravity::Calculate(Tessellation const& tess,
				vector<Primitive> const& cell,int point)
{
  if(_counter==0)
    {
      int n=tess.GetPointNo();
      _length=n;
      masses.resize(n);
      h.resize(n);
      for(int i=0;i<n;++i)
	{
	  masses[i]=cell[i].Density*tess.GetVolume(i);
	  h[i]=tess.GetWidth(i)*_Factor;
	}
      // Create the tree
      ANNkd_tree *tree;
      ANNpointArray treePoints;
      Vector2D pvec;

      treePoints=annAllocPts(n,2);
      for(int i=0;i<n;++i)
	{
	  pvec=tess.GetCellCM(i);
	  treePoints[i][0]=pvec.x;
	  treePoints[i][1]=pvec.y;
	}
      _treePoints=treePoints;
      tree=new ANNkd_tree(masses,_treePoints,n,2,1,ANN_KD_SUGGEST);
      tree->SetH(h);
      _tree=tree;
    }
  Vector2D pvec;
  pvec=tess.GetCellCM(point);
  vector<double> qpoint;
  qpoint.resize(2);
  qpoint[0]=pvec.x;
  qpoint[1]=pvec.y;
  qpoint=_tree->GetAcc(qpoint,_OpenAngle);
  pvec.x=qpoint[0];
  pvec.y=qpoint[1];
  ++_counter;
  if(_counter==_length)
    {
      _counter=0;
      annDeallocPts(_treePoints);
      delete _tree;
      _tree=0;
      annClose();
    }
  return pvec;
}

Vector2D SelfGravity::Calculate
(Tessellation const& tess,
 vector<Primitive> const& cells,
 int point,
 vector<Conserved> const& /*fluxes*/,
 vector<Vector2D> const& /*point_velocity*/,
 HydroBoundaryConditions const& /*hbc*/,
 vector<vector<double> > const& /*tracers*/,
 double /*t*/,
 double /*dt*/)
{
  return Calculate(tess,cells,point);
}