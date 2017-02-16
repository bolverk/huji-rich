#ifndef HDSIM_3D_HPP
#define HDSIM_3D_HPP 1

#include "computational_cell.hpp"
#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "conserved_3d.hpp"
#include "../common/equation_of_state.hpp"
#include "point_motion_3d.hpp"
#include "time_step_function3D.hpp"
#include "flux_calculator_3d.hpp"
#include "cell_updater_3d.hpp"
#include "extensive_updater3d.hpp"

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

    double time;
    size_t cycle;
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
#ifdef RICH_MPI
	  Tessellation3D& tproc,
#endif//RICH_MPI
	  const vector<ComputationalCell3D>& cells,
	  const EquationOfState& eos,
	  const PointMotion3D& pm,
	  const TimeStepFunction3D& tsc,
	  const FluxCalculator3D& fc,
	  const CellUpdater3D& cu,
	  const ExtensiveUpdater3D & eu,
	  const TracerStickerNames tsn);

  //! \brief Advances the simulation in time (first order)
  void timeAdvance();

  void timeAdvance2();

  /*! \brief Access to tessellation
    \return Tessellation
   */
  const Tessellation3D& getTesselation(void) const;

  /*! \brief Access to processor tessellation
  \return Tessellation
  */
#ifdef RICH_MPI
  const Tessellation3D& getProcTesselation(void) const;
#endif

  /*! \brief Access to computational cells
    \return Computational cells
   */
  const vector<ComputationalCell3D>& getCells(void) const;

  const vector<Conserved3D>& getExtensives(void) const;

  Tessellation3D& getTesselation(void);

  /*! \brief Access to computational cells
  \return Computational cells
  */
  vector<ComputationalCell3D>& getCells(void);

  vector<Conserved3D>& getExtensives(void);


  double GetTime(void)const;

  TracerStickerNames GetTracerStickerNames(void)const;

  size_t GetCycle(void)const;

  void SetCycle(size_t cycle);

  void SetTime(double t);

private:
  Tessellation3D& tess_;
#ifdef RICH_MPI
  Tessellation3D& tproc_;
#endif
  const EquationOfState& eos_;
  vector<ComputationalCell3D> cells_;
  vector<Conserved3D> extensive_;
  const PointMotion3D& pm_;
  const TimeStepFunction3D& tsc_;
  const FluxCalculator3D& fc_;
  const CellUpdater3D& cu_;
  const ExtensiveUpdater3D& eu_;
  const TracerStickerNames tsn_;
  ProgressTracker pt_;
};

#endif // HDSIM_3D_HPP
