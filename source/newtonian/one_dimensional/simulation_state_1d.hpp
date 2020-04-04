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

class BoolSpatialDistribution
{
    public:
    
        virtual bool operator()(double x) const = 0;
        
        virtual ~BoolSpatialDistribution(void);
};

class SimulationState1D
{
public:
    
  SimulationState1D
  (const vector<double>& vertices,
   const SpatialDistribution1D& density,
   const SpatialDistribution1D& pressure,
   const SpatialDistribution1D& para_velocity,
   const SpatialDistribution1D& perp_velocity,
   const vector<pair<string, const SpatialDistribution1D*> >& tracers,
   const vector<pair<string, const BoolSpatialDistribution*> >& stickers);

  const vector<double>& getVertices(void) const;

  const vector<ComputationalCell>& getCells(void) const;

  const TracerStickerNames& getTracerStickerNames(void) const;

  void updateVertices(const vector<double>& vertices);

  void updateCells(const vector<ComputationalCell>& cells);
    
private:
  vector<double> vertices_;
  vector<ComputationalCell> cells_;
  TracerStickerNames tsn_;
};

#endif // SIMULATION_STATE_1D_HPP
