#ifndef HDSIM_3D_HPP
#define HDSIM_3D_HPP 1

#include "computational_cell.hpp"
#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "conserved_3d.hpp"
#include "../common/equation_of_state.hpp"
#include "point_motion_3d.hpp"
#include "time_step_calculator.hpp"
#include "flux_calculator_3d.hpp"
#include "cell_updater.hpp"

//! \brief Three dimensional simulation
class HDSim3D
{
public:

  //! \brief Tracks the progress of a simulation
  class ProgressTracker
  {
  public:

    //! \brief Class constructor
    ProgressTracker(void);

    /*! \brief Update the progress tracker
      \param dt Time step
     */
    void update(double dt);

    /*! \brief Returns the current time of the simulation
      \return Time of the simulation
     */
    double getTime(void) const;

    /*! \brief Returns the number of times time advance was called
      \return Cycle number
     */
    double getCycle(void) const;

  private:
    double time_;
    int cycle_;
  };

  /*! \brief Class constructor
    \param tess Tessellation
    \param cells Initial computational cells
    \param eos Equation of state
    \param pm Point motion scheme
    \param tsc Time step calculator
    \param fc Flux calculator
    \param cu Cell updater
   */
  HDSim3D(Tessellation3D& tess,
	  const vector<ComputationalCell>& cells,
	  const EquationOfState& eos,
	  const PointMotion3D& pm,
	  const TimeStepCalculator& tsc,
	  const FluxCalculator& fc,
	  const CellUpdater& cu);

  //! \brief Advances the simulation in time (first order)
  void timeAdvance();

  /*! \brief Access to tessellation
    \return Tessellation
   */
  const Tessellation3D& getTesselation(void) const;

  /*! \brief Access to computational cells
    \return Computational cells
   */
  const vector<ComputationalCell>& getCells(void) const;

private:
  Tessellation3D& tess_;
  const EquationOfState& eos_;
  vector<ComputationalCell> cells_;
  vector<Conserved3D> extensive_;
  const PointMotion3D& pm_;
  const TimeStepCalculator& tsc_;
  const FluxCalculator& fc_;
  const CellUpdater& cu_;
  ProgressTracker pt_;
};

#endif // HDSIM_3D_HPP
