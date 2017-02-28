#ifndef MONOPOLE3D_HPP
#define MONOPOLE3D_HPP 1

#include "ConservativeForce3D.hpp"

class MonopoleSelfGravity3D : public Acceleration3D
{
private:
	const size_t resolution_;
	const double smoothlength_;
public:
	MonopoleSelfGravity3D(size_t resolution,double smoothlength);

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Conserved3D>& fluxes, const double time, TracerStickerNames const& tracerstickernames,
		vector<Vector3D> &acc) const;
};

#endif //MONOPOLE3D_HPP
