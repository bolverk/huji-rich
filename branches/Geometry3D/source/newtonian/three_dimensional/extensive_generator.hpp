#ifndef EXTENSIVE_GENERATOR_HPP
#define EXTENSIVE_GENERATOR_HPP 1

#include "conserved_3d.hpp"
#include "computational_cell.hpp"
#include "../common/equation_of_state.hpp"
#include "../../misc/utils.hpp"
#include "../../3D/GeometryCommon/Tessellation3D.hpp"

Conserved3D calc_intensive(const ComputationalCell& cell,
			 const EquationOfState& eos);

class ExtensiveGenerator: public Index2Member<Conserved3D>
{
public:

  ExtensiveGenerator(const vector<ComputationalCell>& cells,
		     const Tessellation3D& tess,
		     const EquationOfState& eos);

  size_t getLength(void) const;

  Conserved3D operator()(size_t i) const;

private:
  const vector<ComputationalCell>& cells_;
  const Tessellation3D& tess_;
  const EquationOfState& eos_;
};

#endif // EXTENSIVE_GENERATOR_HPP
