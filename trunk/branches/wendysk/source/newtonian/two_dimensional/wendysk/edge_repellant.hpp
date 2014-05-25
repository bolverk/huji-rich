#include "../point_motion.hpp"

class EdgeRepellant: public PointMotion
{
public:

  EdgeRepellant(PointMotion& naive,
		double inner_radius,
		double outer_radius,
		int total_specials);

  Vector2D CalcVelocity(int index, 
			Tessellation const* tess,
			vector<Primitive> const& cells,
			double time);

  vector<Vector2D> calcAllVelocities(Tessellation const* tess,
				     vector<Primitive> const& cells,
				     double time);

private:
  PointMotion& naive_;
  const double inner_radius_;
  const double outer_radius_;
  const int total_specials_;
};
