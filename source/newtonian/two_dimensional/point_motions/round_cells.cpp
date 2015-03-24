#include "round_cells.hpp"

RoundCells::RoundCells
(const PointMotion& pm,
 const EquationOfState& eos,
 double chi,
 double eta):
  pm_(pm), eos_(eos), chi_(chi), eta_(eta) {}

Vector2D RoundCells::calc_dw
(size_t i,
 const Tessellation& tess,
 const vector<ComputationalCell>& cells) const
{
  const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
  const Vector2D s = tess.GetCellCM(static_cast<int>(i));
  const double d = abs(s-r);
  const double R = tess.GetWidth(static_cast<int>(i));
  if(d<0.9*eta_*R)
    return Vector2D(0,0);
  const double c = eos_.dp2c(cells[i].density, cells[i].pressure);
  return chi_*c*(s-r)/d*(d > 1.1*eta_*R ? 1 : 
			 (d-0.9*eta_*R)/(0.2*eta_*R));
}

vector<Vector2D> RoundCells::operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   double time) const
{
  vector<Vector2D> res = pm_(tess,cells,time);
  for(size_t i=0;i<res.size();++i){
    res[i] += calc_dw(i,tess,cells);
  }
  return res;
}
