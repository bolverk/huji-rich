#ifndef CYLINDRICAL_COMPLEMENTARY_HPP
#define CYLINDRICAL_COMPLEMENTARY_HPP 1

#include "../SourceTerm.hpp"
#include "../../../tessellation/geometry.hpp"
#include "../physical_geometry.hpp"

//! \brief Adds necessary correction to cylindrical geometry
class CylindricalComplementary : public SourceTerm
{

public:

	/*! \brief Class constructor
	  \param axis Rotation axis
	 */
	explicit CylindricalComplementary(const Axis& axis);

	vector<Extensive> operator()
		(const Tessellation& tess,
			const PhysicalGeometry& pg_,
			const CacheData& cd,
			const vector<ComputationalCell>& cells,
			const vector<Extensive>& fluxes,
			const vector<Vector2D>& point_velocities,
			const double t,
			TracerStickerNames const& tracerstickernames) const;

private:
	const Axis axis_;
};


#endif // CYLINDRICAL_COMPLEMENTARY_HPP
