#include "simple_cell_updater.hpp"

namespace {
  class CellEdgesGetter: public Index2Member<Edge>
  {
  public:

    CellEdgesGetter(const Tessellation& tess, int n):
      tess_(tess), edge_indices_(tess.GetCellEdges(n)) {}

    size_t getLength(void) const
    {
      return edge_indices_.size();
    }

    Edge operator()(size_t i) const
    {
      return tess_.GetEdge(edge_indices_[i]);
    }

  private:
    const Tessellation& tess_;
    const vector<int> edge_indices_;
  };
}

SimpleCellUpdater::SimpleCellUpdater(void) {}

vector<ComputationalCell> SimpleCellUpdater::operator()
(const Tessellation& tess,
 const PhysicalGeometry& pg,
 const EquationOfState& eos,
 const vector<Extensive>& extensives,
 const vector<ComputationalCell>& old) const
{
  vector<ComputationalCell> res = old;
  for(size_t i=0;i<extensives.size();++i){
    const double volume = pg.calcVolume
      (serial_generate(CellEdgesGetter(tess,static_cast<int>(i))));
    res[i].density = extensives[i].mass/volume;
    res[i].velocity = extensives[i].momentum / extensives[i].mass;
    const double energy = extensives[i].energy/extensives[i].mass - 
      0.5*ScalarProd(res[i].velocity, res[i].velocity);
    res[i].pressure = eos.de2p(res[i].density, energy);
    for(std::map<std::string,double>::const_iterator it =
	  extensives.front().tracers.begin();
	it!=extensives.front().tracers.end();++it)
      res[i].tracers[it->first] = it->second/extensives[i].mass;
  }
  return res;
}
