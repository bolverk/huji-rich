#include "simulation_state_1d.hpp"

SimulationState1D::SimulationState1D
(const vector<double>& vertices,
 const SpatialDistribution1D& density,
 const SpatialDistribution1D& pressure,
 const SpatialDistribution1D& para_velocity,
 const SpatialDistribution1D& perp_velocity,
 const vector<pair<string, const SpatialDistribution1D*> >& tracers,
 const vector<pair<string, const BoolSpatialDistribution*> >& stickers):
  vertices_(vertices),
  cells_(vertices.size()-1),
  tsn_()
{
  for(size_t i=0;i<vertices_.size()-1;++i){
    const double x = 0.5*(vertices_.at(i) + vertices_.at(i+1));
    ComputationalCell& cell = cells_.at(i);
    cell.density = density(x);
    cell.pressure = pressure(x);
    cell.velocity = Vector2D(para_velocity(x),
			      perp_velocity(x));
    for(size_t j=0;j<tracers.size();++j){
      cell.tracers.push_back((*tracers.at(j).second)(x));
      tsn_.tracer_names.push_back(tracers.at(j).first);
    }
    for(size_t j=0;j<stickers.size();++j){
      cell.stickers.push_back((*stickers.at(j).second)(x));
      tsn_.sticker_names.push_back(stickers.at(j).first);
    }
  }
}

const vector<double>& SimulationState1D::getVertices(void) const
{
  return vertices_;
}

const vector<ComputationalCell>& SimulationState1D::getCells(void) const
{
  return cells_;
}

const TracerStickerNames& SimulationState1D::getTracerStickerNames(void) const
{
  return tsn_;
}

void SimulationState1D::updateVertices(const vector<double>& vertices)
{
  vertices_ = vertices;
}

void SimulationState1D::updateCells
(const vector<ComputationalCell>& cells)
{
  cells_ = cells;
}
