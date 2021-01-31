/*! \file simulation_state_1d.hpp
  \author Almog Yalinewich
  \brief Package for computational domain and hydro cells
*/

#ifndef SIMULATION_STATE_1D_HPP
#define SIMULATION_STATE_1D_HPP

#include <utility>
#include <vector>
#include <string>
#include "spatial_distribution1d.hpp"
#include "../two_dimensional/computational_cell_2d.hpp"

using std::pair;
using std::vector;
using std::string;

//! \brief Spatial distribution of boolean values
class BoolSpatialDistribution
{
public:

  /*! \brief Evaluate function
    \param x Position
    \return True or false at that location
   */
  virtual bool operator()(double x) const = 0;
  
  virtual ~BoolSpatialDistribution(void);
};

//! \brief Package for computational domain and hydro cells
class SimulationState1D
{
public:

  /*! \brief Class constructor
    \param vertices Positions of the vertices
    \param density Density profile
    \param pressure Pressure profile
    \param para_velocity Velocity in the parallel direction
    \param perp_velocity Velocity in the perpendicular direction
    \param tracers List of tracer names and profiles
    \param stickers List of sticker names and profiles
   */
  SimulationState1D
  (const vector<double>& vertices,
   const SpatialDistribution1D& density,
   const SpatialDistribution1D& pressure,
   const SpatialDistribution1D& para_velocity,
   const SpatialDistribution1D& perp_velocity,
   const vector<pair<string, const SpatialDistribution1D*> >& tracers,
   const vector<pair<string, const BoolSpatialDistribution*> >& stickers);

  /*! \brief Access to vertices
    \return Positions of vertices
   */
  const vector<double>& getVertices(void) const;

  /*! \brief Access to hydro cells
    \return Hydro cells
   */
  const vector<ComputationalCell>& getCells(void) const;

  /*! \brief Access to tracer and sticker names
    \return Tracer and sticker names
   */
  const TracerStickerNames& getTracerStickerNames(void) const;

  /*! \brief Updates positions of vertices
    \param vertices New vertices
   */
  void updateVertices(const vector<double>& vertices);

  /*! \brief Updates hydro cellls
    \param cells New cells
   */
  void updateCells(const vector<ComputationalCell>& cells);
  
  //! \brief Positions of the vertices
  vector<double> vertices_;
  //! \brief Computational cells
  vector<ComputationalCell> cells_;
  //! \brief Names of stickers and tracers
  TracerStickerNames tsn_;
};

#endif // SIMULATION_STATE_1D_HPP
