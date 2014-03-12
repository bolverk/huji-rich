/*! \file edge_repellant.hpp
  \brief Keeps mesh generating points from leaving the computational domain
  \author Elad Steinberg
 */

#define _USE_MATH_DEFINES
#include "../point_motion.hpp"
#include <cmath>

//! \brief Keeps mesh generating points from leaving the computational domain
class EdgeRepellant: public PointMotion
{
public:

  /*! \brief Class constructor
    \param naive Uncorrected point motion scheme
    \param inner_radius TBA
    \param outer_radius TBA
    \param total_specials TBA
    \param move_inner TBA
    \param inner1 TBA
    \param inner2 TBA
    \param inner3 TBA
    \param n2start TBA
    \param n3start TBA
    \param inner_r1 TBA
    \param inner_r2 TBA
    \param inner_r3 TBA
    \param mass TBA
    \param t0 TBA
    \todo Finish documentation
   */
  EdgeRepellant(PointMotion& naive,
		double inner_radius,
		double outer_radius,
		int total_specials,
		bool move_inner=false,
		int inner1=0,
		int inner2=0,
		int inner3=0,
		int n2start=0,
		int n3start=0,
		double inner_r1=0,
		double inner_r2=0,
		double inner_r3=0,
		double mass=0,
		double t0=0);

  Vector2D CalcVelocity(int index, 
			Tessellation const& tess,
			vector<Primitive> const& cells,
			double time);

  vector<Vector2D> calcAllVelocities(Tessellation const& tess,
				     vector<Primitive> const& cells,
				     double time);

private:
  PointMotion& naive_;
  const double inner_radius_;
  const double outer_radius_;
  const int total_specials_;
  const bool move_inner_;
  const int inner_n1_;
  const int inner_n2_;
  const int inner_n3_;
  const int n2start_;
  const int n3start_;
  const double inner_r1_;
  const double inner_r2_;
  const double inner_r3_;
  const double mass_;
  const double t0_;
  vector<double> init_angles_;
  bool first_time_;
};
