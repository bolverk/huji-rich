/*! \file sb_edge_repellant.hpp
  \brief Keeps the mesh generating points from leaving the computational domain
  \author Almog Yalinewich
 */ 

#include "../two_dimensional/point_motion.hpp"

//! \brief Auxiliary point motion that keeps the points from leaving the computational domain
class SBEdgeRepellant: public PointMotion
{
public:

  /*! \brief Class constructor
    \param naive Uncorrected point motion scheme
    \param inner_radius1 TBA
    \param center1 TBA
    \param inner_radius2 TBA
    \param center2 TBA
    \param outer_radius TBA
    \param total_specials TBA
    \todo Add documentation
   */
  SBEdgeRepellant(PointMotion& naive,
		  double inner_radius1,
		  Vector2D const& center1,
		  double inner_radius2,
		  Vector2D const& center2,
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

  bool should_be_decelerated(Vector2D const& mp,
			     Vector2D const& mv) const;

  PointMotion& naive_;
  const double inner_radius1_;
  const Vector2D center1_;
  const double inner_radius2_;
  const Vector2D center2_;
  const double outer_radius_;
  const int total_specials_;
};
