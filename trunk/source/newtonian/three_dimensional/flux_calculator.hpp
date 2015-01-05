#ifndef FLUX_CALCULATOR_HPP
#define FLUX_CALCULATOR_HPP 1

#include "conserved_3d.hpp"
#include "../../3D/Tessellation/Tessellation3D.hpp"
#include "computational_cell.hpp"

class FluxCalculator
{
public:

  virtual vector<Conserved3D> operator()
  (const Tessellation3D& tess,
     const vector<ComputationalCell>& cells,
     const vector<Vector3D>& point_velocities) const = 0;

  virtual ~FluxCalculator(void);
};

#endif // FLUX_CALCULATOR_HPP
