/*! \file hold_still.hpp
  \brief Keeps a certain portion of the grid stationary
  \author Almog Yalinewich
 */

#ifndef HOLD_STILL_HPP
#define HOLD_STILL_HPP

#include "../point_motion.hpp"

//! \brief Hold some cells still, while letting other move according to another point motion scheme
class HoldStill: public PointMotion
{
public:

  //! \brief Decides which cells fit a criterion
  class Condition
  {
  public:

    /*! \brief Decides if a cetrain cell fit a the criterion
      \param index Cell index
      \param tess Tessellation
      \param cells List of hydrodynamic cells
      \param time Time
      \return True if the cell fits, false otherwise
     */
    virtual bool operator()(int index,
			    const Tessellation& tess,
			    const vector<Primitive>& cells,
			    double time) const = 0;

    virtual ~Condition(void);
  };

  /*! \brief Class constructor
    \param raw Base point motion scheme
    \param cond Says which points should be held still
   */
  HoldStill(PointMotion& raw, const Condition& cond);

  Vector2D CalcVelocity(int index,
			Tessellation const& tessellation,
			vector<Primitive> const& primitives,double time);

  vector<Vector2D> calcAllVelocities(Tessellation const& tess,
				     vector<Primitive> const& cells,
				     double time,vector<CustomEvolution*> &cevolve,
				     const vector<vector<double> >& tracers);

private:
  PointMotion& raw_;
  const Condition& cond_;
};

#endif
