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
	explicit SeveralSources(vector<SourceTerm*> forces);
	//! \brief Class destructor
	~SeveralSources(void);

	vector<Extensive> operator()
		(const Tessellation& tess,
			const PhysicalGeometry& pg,
			const CacheData& cd,
			const vector<ComputationalCell>& cells,
			const vector<Extensive>& fluxes,
			const vector<Vector2D>& point_velocities,
			const double t,
			TracerStickerNames const& tracerstickernames) const;

private:
	vector<SourceTerm*> sources_;
};

#endif //SEVERALSOURCES_HPP
