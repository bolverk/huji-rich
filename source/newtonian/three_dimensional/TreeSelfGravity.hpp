#ifndef TREESG_HPP
#define TREESG_HPP 1

#include "ConservativeForce3D.hpp"

class TreeSelfGravity : public Acceleration3D
{
private:
	const size_t nx_, ny_, nz_;
public:
	mutable vector<double> potential;

	TreeSelfGravity(size_t nx=1, size_t ny=1, size_t nz=1);

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Conserved3D>& fluxes, const double time, TracerStickerNames const& tracerstickernames,
		vector<Vector3D> &acc) const;
};

#endif //TREESG_HPP

