#ifndef HOLD_STILL_HPP
#define HOLD_STILL_HPP

#include "../point_motion.hpp"

class HoldStill: public PointMotion
{
public:

  class Condition
  {
  public:

    virtual bool operator()(int index,
			    const Tessellation& tess,
			    const vector<Primitive>& cells,
			    double time) const = 0;

    virtual ~Condition(void);
  };

  HoldStill(PointMotion& raw, const Condition& cond);

  Vector2D CalcVelocity(int index, 
			Tessellation const& tessellation,
			vector<Primitive> const& primitives,double time);

  vector<Vector2D> calcAllVelocities(Tessellation const& tess,
				     vector<Primitive> const& cells,
				     double time,vector<CustomEvolution*> &cevolve);

private:
  PointMotion& raw_;
  const Condition& cond_;
};

#endif
