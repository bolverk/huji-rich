/*! \file SeveralSources.hpp
  \brief Class for a combination of external sources
  \author Elad Steinberg
*/

#ifndef SEVERALSOURCES_HPP
#define SEVERALSOURCES_HPP 1

#include "../SourceTerm.hpp"
#include <vector>

//! \brief Class for a combination of external sources
class SeveralSources :public SourceTerm
{
public:
	//! \brief Class constructor
	SeveralSources(vector<SourceTerm*> forces);
	//! \brief Class destructor
	~SeveralSources(void);

	Conserved Calculate(Tessellation const* tess,
		vector<Primitive> const& cells,int point,
		vector<Conserved> const& fluxes,vector<Vector2D> const& point_velocity,
		HydroBoundaryConditions const*hbc,
		vector<vector<double> > const &tracer_extensive,vector<double> &dtracer,
		double time,double dt);

private:
	vector<SourceTerm*> sources_;
};

#endif //SEVERALSOURCES_HPP
