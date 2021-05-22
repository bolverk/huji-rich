#ifndef QPOLE3D_HPP
#define QPOLE3D_HPP 1

#include "ConservativeForce3D.hpp"

//! \brief Quadrupolar gravity force term
class QuadrupoleGravity3D : public Acceleration3D
{
private:
	mutable vector<double> edges_;
	const double smoothlength_;
	bool output_;
public:
  /*! \brief List of values of the potential
   */
	mutable vector<double> potential;

  /*! \brief Class constructor
    \param res Resolution
    \param smoothlength Smoothing length
    \param output Logging flag
   */
	QuadrupoleGravity3D(size_t res, double smoothlength,bool output=false);

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Conserved3D>& fluxes, const double time, TracerStickerNames const& tracerstickernames,
		vector<Vector3D> &acc) const override;
};

#endif //QPOLE3D_HPP
