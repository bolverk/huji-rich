#ifndef DIFFUSION_FORCE_HPP
#define DIFFUSION_FORCE_HPP 1
#include "../newtonian/three_dimensional/SourceTerm3D.hpp"
#include "Diffusion.hpp"

class DiffusionForce : public SourceTerm3D
{
public:

    DiffusionForce(Diffusion const& diffusion, std::string const& key): diffusion_(diffusion), key_(key), next_dt_(1e-6 * std::numeric_limits<double>::max()){}

    void operator()(const Tessellation3D& tess,const vector<ComputationalCell3D>& cells,
		const vector<Conserved3D>& fluxes,const vector<Vector3D>& point_velocities, const double t,double dt,
		vector<Conserved3D> &extensives) const;

    double SuggestInverseTimeStep(void)const;
private:
    Diffusion const& diffusion_;
    std::string const key_;
    mutable double next_dt_;
};

#endif