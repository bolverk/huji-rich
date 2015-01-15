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

Conserved CylindricalComplementary::Calculate
(Tessellation const& tess,
 const PhysicalGeometry& pg,
 vector<Primitive> const& cells,
 int point,
 vector<Conserved> const& /*fluxes*/,
 vector<Vector2D> const& /*point_velocity*/,
 HydroBoundaryConditions const& /*hbc*/,
 vector<vector<double> > const& /*tracer*/,
 vector<double>& /*dtracer*/,
 vector<double> const& /*lengthes*/,
 double /*t*/,
 double /*dt*/)
{
  const double p = cells[static_cast<size_t>(point)].Pressure;
  const double r = distance_from_axis(tess.GetCellCM(point),
				      axis_);
  const Vector2D r_hat = cross_z(axis_.direction);
  const double volume = pg.calcVolume
    (serial_generate(CellEdgesGetter(tess,point)));
  return volume*Conserved
    (0,
     (p/r)*r_hat/abs(r_hat),
     0);
}
