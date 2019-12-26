#ifndef LEUPDATE3D_HPP
#define LEUPDATE3D_HPP 1
#include "ConditionExtensiveUpdater3D.hpp"
#include "ConditionActionFlux1.hpp"

class LagrangianExtensiveUpdater3D : public ExtensiveUpdater3D
{
public:
	LagrangianExtensiveUpdater3D(LagrangianFlux3D const& lflux, EquationOfState const& eos, Ghost3D
		const& ghost, const vector<pair<const ConditionExtensiveUpdater3D::Condition3D*, 
		const ConditionExtensiveUpdater3D::Action3D*> >& sequence);

	void operator()(const vector<Conserved3D>& fluxes, const Tessellation3D& tess,
		const double dt, const vector<ComputationalCell3D>& cells, vector<Conserved3D>& extensives, double time,
		TracerStickerNames const& tracerstickernames, const vector<Vector3D>& edge_velocities,
		std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > const& interp_values) const;
private:
	LagrangianFlux3D const& lflux_;
	EquationOfState const& eos_;
	Ghost3D const& ghost_;
	vector<pair<const ConditionExtensiveUpdater3D::Condition3D*, const ConditionExtensiveUpdater3D::Action3D*> >
		const& sequence_;
};
#endif //LEUPDATE_HPP3D