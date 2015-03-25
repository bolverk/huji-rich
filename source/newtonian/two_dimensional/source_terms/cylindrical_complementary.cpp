#include "cylindrical_complementary.hpp"

namespace {
  double distance_from_axis(const Vector2D& point,
			    const Axis& axis)
  {
    const double hypotenuse = abs(point-axis.origin);
    const double side = std::abs(Projection(point-axis.origin,
					    axis.direction));
    return sqrt(pow(hypotenuse,2)-pow(side,2));
  }

  Vector2D cross_z(const Vector2D& v)
  {
    return Vector2D(v.y,-v.x);
  }

  class CellEdgesGetter: public Index2Member<Edge>
  {
  public:
    
    CellEdgesGetter(const Tessellation& tess, int n):
      tess_(tess), edge_indices_(tess.GetCellEdges(n)) {}

    size_t getLength(void) const
    {
      return edge_indices_.size();
    }

    Edge operator()(size_t i) const
    {
      return tess_.GetEdge(edge_indices_[i]);
    }

  private:
    const Tessellation& tess_;
    const vector<int> edge_indices_;
  };
}

CylindricalComplementary::CylindricalComplementary(const Axis& axis):
  axis_(axis) {}

vector<Extensive> CylindricalComplementary::operator()
(const Tessellation& tess,
 const PhysicalGeometry& pg,
 const vector<ComputationalCell>& cells,
 const vector<Extensive>& /*fluxes*/,
 const vector<Vector2D>& /*point_velocities*/,
 const double /*time*/) const
{
  vector<Extensive> res(static_cast<size_t>(tess.GetPointNo()));
  const Vector2D r_hat = cross_z(axis_.direction);
  for(size_t i=0;i<res.size();++i){
    const double p = cells[i].pressure; 
    const double r = distance_from_axis
      (tess.GetCellCM(static_cast<int>(i)),axis_);
    const double volume = pg.calcVolume
      (serial_generate(CellEdgesGetter(tess,static_cast<int>(i))));
    res[i].mass = 0;
    res[i].momentum = volume*(p/r)*r_hat/abs(r_hat);
    res[i].energy = 0;
  }
  return res;
}
