#include "SourceTerm3D.hpp"
#include <limits>

SourceTerm3D::~SourceTerm3D(void){}

double SourceTerm3D::SuggestInverseTimeStep(void)const
{
	return 100 * std::numeric_limits<double>::min();
}

void ZeroForce3D::operator()(const Tessellation3D& /*tess*/, const vector<ComputationalCell3D>& /*cells*/,
		const vector<Conserved3D>& /*fluxes*/, const vector<Vector3D>& /*point_velocities*/, const double /*t*/, 
		double /*dt*/, vector<Conserved3D> &/*extensives*/) const {} 
