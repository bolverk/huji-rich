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
#include "../../mpi/ProcessorUpdate3D.hpp"
#include "SourceTerm3D.hpp"

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
    void updateTime(double dt);

	void updateCycle();

    /*! \brief Returns the current time of the simulation
      \return Time of the simulation
     */
    double getTime(void) const;

    /*! \brief Returns the number of times time advance was called
      \return Cycle number
     */
    size_t getCycle(void) const;

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
	\param eu Extensive updater
	\param source Source term
	\param tsn The names of the stickers and tracers
	\param proc_update How to load balance
	\param tproc The tessellation of the domian decomposition
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
	  const ExtensiveUpdater3D& eu,
	  const	SourceTerm3D& source,
	  TracerStickerNames& tsn,
	  bool SR=false
#ifdef RICH_MPI
	  ,const ProcessorUpdate3D* proc_update = 0
#endif
	  ,bool new_start = true
	  ,const double maxload = 4.0
  );

  //! \brief Advances the simulation in time (first order)
  void timeAdvance();
  //! \brief Advances the simulation in time (second order)
  void timeAdvance2();

  void timeAdvance3();

  void timeAdvance32();

  void timeAdvance33();

  void timeAdvance4();

  /*! \brief Access to tessellation
    \return Tessellation
   */
  const Tessellation3D& getTesselation(void) const;

  /*! \brief Access to processor tessellation
  \return Tessellation
  */
#ifdef RICH_MPI
  const Tessellation3D& getProcTesselation(void) const;

  Tessellation3D& getProcTesselation(void);
#endif

  /*! \brief Access to computational cells
    \return Computational cells
   */
  const vector<ComputationalCell3D>& getCells(void) const;
  /*! \brief Access to extensive cells
  \return Extensive cells
  */
  const vector<Conserved3D>& getExtensives(void) const;
  /*! \brief Access to tessellation
  \return Tessellation
  */
  Tessellation3D& getTesselation(void);

  /*! \brief Access to computational cells
  \return Computational cells
  */
  vector<ComputationalCell3D>& getCells(void);
  /*! \brief Access to extensive cells
  \return Extensive cells
  */
  vector<Conserved3D>& getExtensives(void);
  /*! \brief Get the time of the simulation
  \return The current time of the simulation
  */
  double GetTime(void)const;
  /*! \brief Access to the names of the stickers and tracers
  \return The names of the stickers and tracers
  */
  TracerStickerNames GetTracerStickerNames(void)const;
  /*! \brief Get the cycle number of the simulation
  \return The current cycle of the simulation
  */
  size_t GetCycle(void)const;
  /*! \brief Change the cycle of the simulation
  \param cycle The new cycle of the simulation
  */
  void SetCycle(size_t cycle);
  /*! \brief Change the time of the simulation
  \param t The new time of the simulation
  */
  void SetTime(double t);

  /*! \brief Change/get the max ID of the sim
  \return t The maximum ID number ofr all cells
  */
  size_t & GetMaxID(void);

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
  const	SourceTerm3D &source_;
  TracerStickerNames &tsn_;
  ProgressTracker pt_;
#ifdef RICH_MPI
  const ProcessorUpdate3D* proc_update_;
#endif
  size_t Max_ID_;
  const double maxload_;
};

#endif // HDSIM_3D_HPP
