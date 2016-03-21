/*! \file ConservativeForce.hpp
  \brief Abstract class for conservative force's acceleration
  \author Elad Steinberg
*/

#ifndef CONSFORCE_HPP
#define CONSFORCE_HPP 1

#include "../SourceTerm.hpp"

//! \brief Physical acceleration
class Acceleration
{
public:
	/*!
	  \brief Calculates the acceleration that the cell feels
	  \param tess The tessellation
	  \param cells The primitive cells
	  \param point The index of the cell to calculate
	  \param fluxes The vector of the fluxes
	  \param time The simulation time
	  \param tracerstickernames The names of the tracers and stickers
	  \return The calculated acceleration
	*/
	virtual Vector2D operator()
		(const Tessellation& tess,
			const vector<ComputationalCell>& cells,
			const vector<Extensive>& fluxes,
			const double time,
			const int point,
			TracerStickerNames const& tracerstickernames) const = 0;

	virtual ~Acceleration(void);
};
/*! \brief Class for conservative forces
  \author Elad Steinberg
*/
class ConservativeForce : public SourceTerm
{
public:
	/*! \brief Class constructor
	  \param acc The acceleration force
	  \param mass_flux Flag whether to include mass flux into energy equation
	*/
	explicit ConservativeForce(const Acceleration& acc, bool mass_flux = false);

	/*!
	  \brief Class destructor
	*/
	~ConservativeForce(void);

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
	const Acceleration& acc_;
	const bool mass_flux_;

	ConservativeForce(const ConservativeForce& origin);
	ConservativeForce& operator=(const ConservativeForce& origin);
};

#endif // CONSFORCE_HPP
