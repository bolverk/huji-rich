#ifndef ANNTREESG_HPP
#define ANNTREESG_HPP 1

#include "ConservativeForce3D.hpp"

class ANNSelfGravity : public Acceleration3D
{
private:
	Tessellation3D const* tproc_;
	const double opening_;
public:
	ANNSelfGravity(double opening = 0.25,Tessellation3D const* tproc=0);

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Conserved3D>& fluxes, const double time, TracerStickerNames const& tracerstickernames,
		vector<Vector3D> &acc) const;
};

#endif //ANNTREESG_HPP

