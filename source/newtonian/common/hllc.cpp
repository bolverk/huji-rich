#include "hllc.hpp"
#include "hydrodynamic_variables.hpp"
#include "../../misc/universal_error.hpp"
#include "../../misc/utils.hpp"

using namespace std;

namespace {
  class WaveSpeeds
  {
  public:

    WaveSpeeds(double left_i,
	       double center_i,
	       double right_i):
      left(left_i),
      center(center_i),
      right(right_i) {}

    const double left;
    const double center;
    const double right;
  };
}

namespace {
WaveSpeeds estimate_wave_speeds
(Primitive const& left, Primitive const& right)
{
  const double dl = left.Density;
  const double pl = left.Pressure;
  const double vl = left.Velocity.x;
  const double cl = left.SoundSpeed;
  const double dr = right.Density;
  const double pr = right.Pressure;
  const double vr = right.Velocity.x;
  const double cr = right.SoundSpeed;
  const double sl = min(vl-cl,vr-cr);
  const double sr = max(vl+cl,vr+cr);
  const double ss = (pr-pl+dl*vl*(sl-vl)-dr*vr*(sr-vr))/
    (dl*(sl-vl)-dr*(sr-vr));
  return WaveSpeeds(sl,ss,sr);
}
}

namespace {
  UniversalError invalid_wave_speeds(Primitive const& left,
				     Primitive const& right,
				     double velocity,
				     double left_wave_speed,
				     double center_wave_speed,
				     double right_wave_speed)
  {
    UniversalError res("Invalid wave speeds in hllc solver");
    res.AddEntry("left density",left.Density);
    res.AddEntry("left pressure",left.Pressure);
    res.AddEntry("left x velocity",left.Velocity.x);
    res.AddEntry("left y velocity",left.Velocity.y);
    res.AddEntry("left sound speed",left.SoundSpeed);
    res.AddEntry("left energy",left.Energy);
    res.AddEntry("right density",right.Density);
    res.AddEntry("right pressure",right.Pressure);
    res.AddEntry("right x velocity",right.Velocity.x);
    res.AddEntry("right y velocity",right.Velocity.y);
    res.AddEntry("right sound speed",right.SoundSpeed);
    res.AddEntry("right energy",right.Energy);
    res.AddEntry("interface velocity",velocity);
    res.AddEntry("left wave speed",left_wave_speed);
    res.AddEntry("center wave speed",center_wave_speed);
    res.AddEntry("right wave speed",right_wave_speed);
    return res;
  }
}

namespace {
Conserved starred_state
(Primitive const& state, double sk, double ss)
{
  const double dk = state.Density;
  const double pk = state.Pressure;
  const double uk = state.Velocity.x;
  const double vk = state.Velocity.y;
  const double ds = dk*(sk-uk)/(sk-ss);
  const double ek = TotalEnergyDensity(state);
  Conserved res;
  res.Mass = ds;
  res.Momentum.x = ds*ss;
  res.Momentum.y = ds*vk;
  res.Energy = ek*ds/dk+
    ds*(ss-uk)*(ss+pk/dk/(sk-uk));
  return res;
}
}

Conserved Hllc::operator()
  (Primitive const& left,
   Primitive const& right,
   double velocity) const
{
  if(is_nan(right.Velocity.x))
    throw UniversalError("Hllc::Solved entered with nan");

  const Vector2D normaldir(1,0);
  Primitive local_left=left;
  Primitive local_right=right;

  local_left.Velocity -= velocity*normaldir;
  local_right.Velocity -= velocity*normaldir;

  const Conserved ul = Primitive2Conserved(local_left);
  const Conserved ur = Primitive2Conserved(local_right);

  const Vector2D xdir(1,0);
  const Conserved fl = Primitive2Flux(local_left, xdir);
  const Conserved fr = Primitive2Flux(local_right, xdir);

  const WaveSpeeds ws = estimate_wave_speeds(local_left, local_right);

  const Conserved usl = starred_state(local_left, ws.left, ws.center);
  const Conserved usr = starred_state(local_right, ws.right, ws.center);

  Conserved f_gr;
  if(ws.left>0)
    f_gr = fl;
  else if(ws.left<=0&&ws.center>=0)
   f_gr = fl+ws.left*(usl-ul);
  else if(ws.center<0&&ws.right>=0)
    f_gr  = fr + ws.right*(usr-ur);
  else if(ws.right<0)
    f_gr = fr;
  else
    throw invalid_wave_speeds(local_left,
			      local_right,
			      velocity,
			      ws.left,
			      ws.center,
			      ws.right);

  f_gr.Energy += ScalarProd(f_gr.Momentum,velocity*normaldir) +
                0.5*f_gr.Mass*velocity*velocity;
  f_gr.Momentum += velocity*f_gr.Mass*normaldir;
  return f_gr;
}
